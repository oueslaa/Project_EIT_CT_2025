function Batch(varargin)
%BATCH  Relance safe & idempotent.
% - Skip si outputs déjà présents (eit_pack/mesh/E_electrodes)
% - Ne refait que le Plot si le Forward est déjà fait mais pas les PNG
% - Checkpoint dans BatchLogs/completed.tsv

%% --- Args (facultatifs) -------------------------------------------------
p = inputParser;
addParameter(p,'RequirePlots',false,@islogical);            % si true: "done" seulement si PNGs clés existent
addParameter(p,'OnlyPatients',{},@(x) iscellstr(x) || ischar(x));
addParameter(p,'OnlyZSlices',[],@(x) isempty(x) || isnumeric(x));  % ex: 280:300 pour debug
parse(p,varargin{:});
RequirePlots = p.Results.RequirePlots;
if ischar(p.Results.OnlyPatients), onlyPatients = {p.Results.OnlyPatients}; else, onlyPatients = p.Results.OnlyPatients; end
onlyZSlices  = p.Results.OnlyZSlices;

%% --- Setup --------------------------------------------------------------
projRoot = fileparts(mfilename('fullpath')); cd(projRoot); addpath(projRoot);
addpath('src_anis', genpath('src_anis'));

ooeitDir = fullfile(fileparts(projRoot), 'OOEIT-main');
if exist(ooeitDir,'dir'), addpath(genpath(ooeitDir)); end

% Logs / checkpoint
logdir   = fullfile(projRoot, 'BatchLogs'); if ~exist(logdir,'dir'), mkdir(logdir); end
ts       = datestr(now,'yyyymmdd_HHMMSS');
logfile  = fullfile(logdir, [ts '_BatchTorsoForward.log']);
ckptFile = fullfile(logdir, 'completed.tsv');   % (patient \t z)

% Charger le checkpoint existant (si présent) – lecture ligne à ligne robuste
doneSet = containers.Map('KeyType','char','ValueType','logical');
if exist(ckptFile,'file')==2
    fid = fopen(ckptFile,'r');
    if fid>0
        tline = fgetl(fid);
        while ischar(tline)
            if ~isempty(tline)
                parts = strsplit(strtrim(tline), sprintf('\t'));
                if numel(parts) >= 2
                    pid = parts{1};
                    z   = str2double(parts{2});
                    if ~isnan(z)
                        key = sprintf('%s:%d', pid, z);
                        doneSet(key) = true;
                    end
                end
            end
            tline = fgetl(fid);
        end
        fclose(fid);
    end
end


diary(logfile); diary on;
fprintf('== Batch Torso Forward (start %s) ==\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));

%% --- Lire meta.csv & filtrer --------------------------------------------
metaFile = fullfile(projRoot, 'meta', 'meta.csv');
assert(exist(metaFile,'file')==2, 'meta/meta.csv introuvable');

opts = detectImportOptions(metaFile, 'Delimiter',';', 'TextType','char');  % pas de string
T = readtable(metaFile, opts);

colImage = find(strcmpi(T.Properties.VariableNames,'image_id'), 1);
colStudy = find(strcmpi(T.Properties.VariableNames,'study_type'), 1);
assert(~isempty(colImage) && ~isempty(colStudy), 'Colonnes image_id/study_type manquantes');

image_ids = T{:, colImage};         % cellstr ou char -> on force cellstr après
if ischar(image_ids), image_ids = cellstr(image_ids); end
study_vec = T{:, colStudy};         % idem
if ischar(study_vec), study_vec = cellstr(study_vec); end

mask = false(numel(study_vec),1);
for i=1:numel(study_vec)
    mask(i) = strcmpi(strtrim(study_vec{i}), 'ct neck-thorax-abdomen-pelvis');
end
sel_ids = image_ids(mask);

% filtre optionnel de patients
if ~isempty(onlyPatients)
    keep = false(numel(sel_ids),1);
    for i=1:numel(sel_ids)
        keep(i) = any(strcmp(sel_ids{i}, onlyPatients));
    end
    sel_ids = sel_ids(keep);
end

fprintf('[Batch] %d sujets gardés sur %d (study_type=neck-thorax-abdomen-pelvis)\n', numel(sel_ids), height(T));

% Patterns des fichiers de segmentation
heart_pat = '(heart|atrial[_-]?appendage)';
lung_pat  = '(lung|pulmonary[_-]?lobe|upper[_-]?lobe|middle[_-]?lobe|lower[_-]?lobe)';

%% --- Boucle sujets -------------------------------------------------------
for ii = 1:numel(sel_ids)
    patient_id = sel_ids{ii};
    try
        data_dir = fullfile(projRoot, 'Data_set', patient_id);
        seg_dir  = fullfile(data_dir, 'segmentations');
        if ~isfolder(seg_dir)
            warning('[%s] segmentations/ introuvable -> skip', patient_id);
            continue;
        end

        zs = find_slices_heart_lungs(seg_dir, heart_pat, lung_pat);
        if ~isempty(onlyZSlices), zs = intersect(zs, onlyZSlices); end

        if isempty(zs)
            fprintf('[%s] Aucune slice (coeur & poumons) -> skip\n', patient_id);
            continue;
        end
        fprintf('[%s] %d slice(s): z={%s}\n', patient_id, numel(zs), numlist(zs));

        for z_slice = zs(:).'
            key = sprintf('%s:%d', patient_id, z_slice);

            % checkpoint + état disque
            if isKey(doneSet, key) && slice_is_done(patient_id, z_slice, RequirePlots)
                fprintf('[%s][z=%d] déjà fait (checkpoint) -> skip\n', patient_id, z_slice);
                continue;
            end

            [fwdDone, plotsDone] = slice_state(patient_id, z_slice);

            try
                if ~fwdDone
                    Problem_Forward_fn(patient_id, z_slice);
                    fwdDone = true;
                end

                if ~plotsDone
                    Plot_Problem_Forward_fn(patient_id, z_slice, ...
                        'SavePNG', true, 'DPI', 300, 'ShowFigures', false);
                    plotsDone = true;
                end

                if fwdDone && (~RequirePlots || plotsDone)
                    append_checkpoint(ckptFile, patient_id, z_slice);
                    doneSet(key) = true;
                    fprintf('[%s][z=%d] OK\n', patient_id, z_slice);
                else
                    warning('[%s][z=%d] incomplet: fwd=%d, plots=%d', patient_id, z_slice, fwdDone, plotsDone);
                end

            catch MEz
                warning('[%s][z=%d] échec: %s', patient_id, z_slice, MEz.message);
            end
        end

    catch MEsub
        warning('[%s] erreur sujet: %s', patient_id, MEsub.message);
    end
end

fprintf('== Batch terminé (%s) ==\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
diary off;
end

%% =========================== Helpers =====================================
function zs = find_slices_heart_lungs(seg_dir, heart_pat, lung_pat)
    S = dir(fullfile(seg_dir, '*.nii*'));
    if isempty(S), zs = []; return; end
    heart_z = []; lung_z = [];
    for i = 1:numel(S)
        nm = lower(S(i).name);
        fp = fullfile(seg_dir, S(i).name);
        try
            if ~isempty(regexp(nm, heart_pat, 'once'))
                V = read_nifti3D_flex(fp); heart_z = union(heart_z, nonzero_slices(V));
            elseif ~isempty(regexp(nm, lung_pat, 'once'))
                V = read_nifti3D_flex(fp); lung_z  = union(lung_z,  nonzero_slices(V));
            end
        catch
            % ignorer masque corrompu
        end
    end
    zs = intersect(heart_z, lung_z); zs = zs(:).';
end

function V = read_nifti3D_flex(fp)
    if exist('read_nifti3D','file')==2
        try, [V,~] = read_nifti3D(fp); return; catch, end
    end
    V = niftiread(fp);
end

function zs = nonzero_slices(V)
    if ndims(V)~=3, zs=[]; return; end
    V = logical(V);
    nz = squeeze(any(any(V,1),2));
    zs = find(nz(:)>0);
end

function tf = slice_is_done(patient_id, z_slice, requirePlots)
    [fwdDone, plotsDone] = slice_state(patient_id, z_slice);
    tf = fwdDone && (~requirePlots || plotsDone);
end

function [fwdDone, plotsDone] = slice_state(patient_id, z_slice)
    rootOut = fullfile('Outputs', patient_id, sprintf('slice_%03d', z_slice));
    meshFile = fullfile(rootOut, 'mesh', sprintf('mesh_slice%03d.mat', z_slice));
    packFile = fullfile(rootOut, 'eit_pack.mat');
    elFile   = fullfile(rootOut, 'E_electrodes.mat');

    fwdDone = exist(meshFile,'file')==2 && exist(packFile,'file')==2 && exist(elFile,'file')==2;

    plotDir  = fullfile(rootOut,'plots_forward');
    keyPNGs = { 'ct_brut.png', ...
                'mesh_electrodes_numbered.png', ...
                'mesh_electrodes_sigma_forward.png', ...
                'electrode_voltages_heatmap.png', ...
                'reciprocity_heatmap.png' };
    if ~exist(plotDir,'dir')
        plotsDone = false;
    else
        ex = false(numel(keyPNGs),1);
        for i=1:numel(keyPNGs)
            ex(i) = exist(fullfile(plotDir,keyPNGs{i}),'file')==2;
        end
        plotsDone = all(ex);
    end
end

function append_checkpoint(ckptFile, patient_id, z_slice)
    fid = fopen(ckptFile,'a');
    if fid<0
        warning('Impossible d''ouvrir le checkpoint: %s', ckptFile);
        return;
    end
    fprintf(fid, '%s\t%d\n', patient_id, z_slice);
    fclose(fid);
end

function s = numlist(zs)
    % zs -> 'a, b, c' (char)
    if isempty(zs), s = ''; return; end
    zs = zs(:).';           % row
    s = sprintf('%d, ', zs);
    s = s(1:end-2);
end
