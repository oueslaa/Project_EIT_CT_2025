% ==========================================
% File: Plot_Problem_Inverse_ConstInit.m
% ==========================================
% Affiche TOUT (réf / init DISCRÈTE / itérations / reco / montages + métriques)
% pour une reco faite avec initialisation homogène (sigma constant).

clear; close all; clc;
addpath(genpath('src_anis'));
addpath(genpath('/Users/anis/Documents/StageInria/Code/OOEIT-main'));  % adapte si besoin
set(groot,'defaultFigureVisible','on');

% --- Sauvegarde PNG ?
SAVE_PNG = true;  DPI = 300;

% --- Paramètres utilisateur ---
patient_id = 's0011';
z_slice    = 301;

% (facultatif) valeur attendue ; la détection auto corrige si besoin
INIT_CONST = 0.35;   % <- si faux, le script retrouvera le bon run
num2tag = @(x) strrep(sprintf('%.4g',x),'.','p');

% --- Dossiers de base ---
rootOut = fullfile('Outputs', patient_id, sprintf('slice_%03d', z_slice));  % pack/mesh
invBase = fullfile(rootOut, 'inverse');
invRoot = fullfile(invBase, sprintf('initConst_%s', num2tag(INIT_CONST)));
recFile = fullfile(invRoot, 'reconstruction_inverse.mat');

% Si le fichier n'existe pas, on cherche automatiquement le dernier run initConst_*
if ~isfile(recFile)
    C = dir(fullfile(invBase,'initConst_*','reconstruction_inverse.mat'));
    if ~isempty(C)
        [~,ix] = max([C.datenum]);
        invRoot = C(ix).folder;
        recFile = fullfile(invRoot,'reconstruction_inverse.mat');
        fprintf('[PLOT-CONST] initConst introuvable pour %.4g, prise du run le plus récent: %s\n', ...
            INIT_CONST, invRoot);
    else
        error('Aucune reconstruction initConst_* trouvée dans %s', invBase);
    end
end

plotDir = fullfile(invRoot, 'plots_inverse'); if ~exist(plotDir,'dir'), mkdir(plotDir); end
iterDir = fullfile(plotDir, 'iterations');

packFile = fullfile(rootOut, 'eit_pack.mat');       % pack forward (facultatif pour GT/viz)
vizFile  = fullfile(rootOut, 'viz_config.mat');

% --- Chargements ---
R = load(recFile);   % sigma0, sigma_rec, sigma_rec_tri, initInfo, g,H, misfit, Umeas, Ufwd, scaleU, ...
S = []; if isfile(packFile), S = load(packFile); end
viz = []; if isfile(vizFile), tmp = load(vizFile); if isfield(tmp,'viz_config'), viz = tmp.viz_config; end, end

% Valeur réelle utilisée (si dispo)
if isfield(R,'initInfo') && isfield(R.initInfo,'levels') && numel(R.initInfo.levels)==1
    INIT_CONST = R.initInfo.levels(1);
end

% géométrie / maillage
if ~isempty(S) && isfield(S,'g') && isfield(S,'H')
    g=S.g; H=S.H;
else
    g=R.g; H=R.H;
end
assert(~isempty(g) && ~isempty(H), 'Maillage introuvable');

% mode d'affichage
displayMode = 'neurological';
if ~isempty(viz) && isfield(viz,'displayMode'), displayMode = viz.displayMode; end

% ---------- Référence ----------
triGT = [];
if isfield(R,'sigma_true_tri') && ~isempty(R.sigma_true_tri)
    triGT = R.sigma_true_tri(:);
elseif ~isempty(S) && isfield(S,'sigma_tri') && ~isempty(S.sigma_tri)
    triGT = S.sigma_tri(:);
end

% ---------- Init DISCRÈTE (sera constante) ----------
if isfield(R,'initInfo') && isfield(R.initInfo,'lab_tri') && isfield(R.initInfo,'levels') ...
        && ~isempty(R.initInfo.lab_tri) && ~isempty(R.initInfo.levels)
    triInit = R.initInfo.levels(R.initInfo.lab_tri(:));
else
    levels_guess = infer_levels_from_packs(R, S);
    triInit = tri_from_nodal_mode(H, R.sigma0, levels_guess);
end
assert(~isempty(triInit), 'Init manquante');

% ---------- Reco ----------
if isfield(R,'sigma_rec_tri') && ~isempty(R.sigma_rec_tri)
    triRec  = R.sigma_rec_tri(:);
else
    triRec  = mean([R.sigma_rec(H(:,1)) R.sigma_rec(H(:,2)) R.sigma_rec(H(:,3))],2);
end

% ---------- Echelle couleur ----------
climUsed = [];
if ~isempty(viz) && isfield(viz,'sigma_clim') && numel(viz.sigma_clim)==2 && all(isfinite(viz.sigma_clim))
    climUsed = viz.sigma_clim;
elseif isfield(R,'sigma_clim') && numel(R.sigma_clim)==2 && all(isfinite(R.sigma_clim))
    climUsed = R.sigma_clim;
else
    vmin = min([triInit; triRec]); vmax = max([triInit; triRec]);
    if ~isempty(triGT), vmin=min(vmin,min(triGT)); vmax=max(vmax,max(triGT)); end
    climUsed = [vmin vmax];
end

%% ---------- Plot 1: Référence ----------
if ~isempty(triGT) && numel(triGT)==size(H,1)
    f1 = figure('Color','w','Name','Sigma - Référence','NumberTitle','off');
    tri_colorplot(g,H,triGT,climUsed,displayMode);
    title(sprintf('Cible (\\sigma_{true}) — %s slice %d [initConst=%.3g]', patient_id, z_slice, INIT_CONST));
    cb=colorbar; ylabel(cb,'S/m');
    if SAVE_PNG, exportgraphics(f1, fullfile(plotDir, sprintf('sigma_ref_slice%03d_initConst_%s.png',z_slice,num2tag(INIT_CONST))), 'Resolution', DPI); end
elseif ~isempty(triGT)
    warning('GT de taille %d incompatible avec Ntri=%d : ignorée.', numel(triGT), size(H,1));
else
    fprintf('Ground truth indisponible : skip plot référence\n');
end

%% ---------- Plot 2: Init (discrète) ----------
f2 = figure('Color','w','Name','Sigma - Init (discrète)','NumberTitle','off');
tri_colorplot(g,H,triInit,climUsed,displayMode);
title(sprintf('Init (\\sigma_0=%.3g) — %s slice %d', INIT_CONST, patient_id, z_slice));
cb=colorbar; ylabel(cb,'S/m');
if SAVE_PNG, exportgraphics(f2, fullfile(plotDir, sprintf('sigma_init_slice%03d_initConst_%s.png',z_slice,num2tag(INIT_CONST))), 'Resolution', DPI); end

%% ---------- Intermédiaires par STAGE (tous les snapshots du run) ----------
if exist(iterDir,'dir')
    D = dir(fullfile(iterDir,'iter_*.mat'));
    if ~isempty(D)
        [~,ord] = sort({D.name}); D = D(ord);
        stage_ids = zeros(1,numel(D));
        for i = 1:numel(D)
            Si = load(fullfile(iterDir, D(i).name));
            if isfield(Si,'stage_k') && ~isempty(Si.stage_k) && isnumeric(Si.stage_k)
                stage_ids(i) = Si.stage_k;
            else
                stage_ids(i) = 1;
            end
        end
        stages = unique(stage_ids);
        for sID = stages
            idx = find(stage_ids==sID); n = numel(idx);
            ncols = max(1, ceil(sqrt(n))); nrows = max(1, ceil(n / ncols));
            fS  = figure('Color','w','Name',sprintf('Reconstructions intermédiaires — stage %d',sID),'NumberTitle','off');
            tlo = tiledlayout(fS, nrows, ncols, 'Padding','compact','TileSpacing','compact');

            axFirst = [];
            for k = 1:n
                ax = nexttile(tlo);
                if isempty(axFirst), axFirst = ax; end
                Si = load(fullfile(iterDir, D(idx(k)).name));
                if isfield(Si,'sigma_iter_tri') && ~isempty(Si.sigma_iter_tri)
                    triIt = Si.sigma_iter_tri(:);
                else
                    if ~isfield(Si,'sigma_iter_nodal') || isempty(Si.sigma_iter_nodal), continue; end
                    nod   = Si.sigma_iter_nodal(:);
                    triIt = mean([nod(H(:,1)) nod(H(:,2)) nod(H(:,3))],2);
                end
                if numel(triIt) ~= size(H,1)
                    warning('Snapshot %s ignoré (taille tri %d vs Ntri %d).', D(idx(k)).name, numel(triIt), size(H,1));
                    continue;
                end
                tri_colorplot(g,H,triIt,climUsed,displayMode);
                ttl = 'iter';
                if isfield(Si,'iter_global') && ~isempty(Si.iter_global), ttl = sprintf('iter %d', Si.iter_global); end
                title(ttl);
            end
            if ~isempty(axFirst)
                cb = colorbar(axFirst);
                try, cb.Layout.Tile = 'east'; catch, set(cb,'Location','eastoutside'); end
                ylabel(cb,'S/m');
            end
            if SAVE_PNG
                exportgraphics(fS, fullfile(plotDir, ...
                    sprintf('intermediate_recos_stage%d_slice%03d_initConst_%s.png',sID,z_slice,num2tag(INIT_CONST))), ...
                    'Resolution', DPI);
            end
        end
    else
        fprintf('Aucun snapshot trouvé dans %s\n', iterDir);
    end
else
    fprintf('Répertoire des itérations absent: %s\n', iterDir);
end

%% ---------- Plot 3: Reco ----------
f3 = figure('Color','w','Name','Sigma - Reco','NumberTitle','off');
tri_colorplot(g,H,triRec,climUsed,displayMode);
mis = NaN; if isfield(R,'misfit'), mis = R.misfit; end
title(sprintf('Reco — misfit=%.4g — %s slice %d [initConst=%.3g]', mis, patient_id, z_slice, INIT_CONST));
cb=colorbar; ylabel(cb,'S/m');
if SAVE_PNG, exportgraphics(f3, fullfile(plotDir, sprintf('sigma_reco_slice%03d_initConst_%s.png',z_slice,num2tag(INIT_CONST))), 'Resolution', DPI); end

%% ---------- Plot 4: Montage 1×3 ----------
f4 = figure('Color','w','Name','Montage 1x3 — Réf / Init / Reco','NumberTitle','off');
tiledlayout(f4,1,3,'Padding','compact','TileSpacing','compact');

% Réf
nexttile;
if ~isempty(triGT) && numel(triGT)==size(H,1)
    tri_colorplot(g,H,triGT,climUsed,displayMode);
    title('Cible (\sigma_{true})'); cb=colorbar; ylabel(cb,'S/m');
else
    text(0.5,0.5,'Ground truth indisponible','Units','normalized', ...
         'HorizontalAlignment','center'); axis off;
end

% Init (discrète)
nexttile;
tri_colorplot(g,H,triInit,climUsed,displayMode);
title(sprintf('Init (\\sigma_0=%.3g)', INIT_CONST)); cb=colorbar; ylabel(cb,'S/m');

% Reco
nexttile;
tri_colorplot(g,H,triRec,climUsed,displayMode);
title(sprintf('Reco (misfit=%.4g)', mis)); cb=colorbar; ylabel(cb,'S/m');

if SAVE_PNG, exportgraphics(f4, fullfile(plotDir, sprintf('montage_ref_init_reco_slice%03d_initConst_%s.png',z_slice,num2tag(INIT_CONST))), 'Resolution', DPI); end

%% ---------- (Optionnel) Métriques & trajectoires ----------
try
    if exist('EITFEM','class') == 8 && isfield(S,'E') && isfield(S,'params') ...
            && isfield(R,'Umeas') && ~isempty(R.Umeas)
        fem = EITFEM(ForwardMesh1st(g*1e-3, H, S.E)); fem.mode='current';
        fem.zeta = S.params.z_contact*ones(S.params.Ne,1);
        Umeas = R.Umeas(:); Ne = S.params.Ne; Kmeas = numel(Umeas)/Ne;

        if isfield(S,'Ipat') && ~isempty(S.Ipat)
            base=S.Ipat;
        else
            if exist('EITSim','class') == 8 || exist('EITSim','class') == 2
                base=EITSim.buildTrigPattern(Ne, S.params.I_amp);
            else
                base=inj_pattern_trig(Ne, S.params.I_amp);
            end
        end
        if size(base,2) < Kmeas, base = repmat(base,1,ceil(Kmeas/size(base,2))); end
        inj = base(:,1:Kmeas); fem.Iel = inj(:); fem.Uel = Umeas; fem.SetInvGamma(5e-6,8e-3);

        D = dir(fullfile(iterDir,'iter_*.mat'));
        [~,ord]=sort({D.name}); D=D(ord);
        vals = nan(numel(D),1);
        for i=1:numel(D)
            Si=load(fullfile(iterDir,D(i).name));
            if isfield(Si,'sigma_iter_nodal') && ~isempty(Si.sigma_iter_nodal)
                sig = Si.sigma_iter_nodal(:);
            elseif isfield(Si,'sigma_iter_tri') && ~isempty(Si.sigma_iter_tri)
                sig = nodal_from_tri(g,H,Si.sigma_iter_tri(:));
            else
                continue;
            end
            U  = fem.SolveForwardVec(sig);
            vals(i) = norm(U - Umeas) / max(norm(Umeas), eps);
        end
        if ~isempty(vals) && any(isfinite(vals))
            fL2 = figure('Color','w','Name','L2(U_{iter}-U_{meas})','NumberTitle','off');
            plot(vals,'-o'); grid on; xlabel('itération (ordre fichiers)'); ylabel('L2 relatif');
            title('Convergence data-fit');
            if SAVE_PNG
                exportgraphics(fL2, fullfile(plotDir, sprintf('l2_datafit_slice%03d_initConst_%s.png',z_slice,num2tag(INIT_CONST))), 'Resolution', DPI);
            end
        end
    end
catch ME
    warning('Plot L2 data-fit sauté: %s', ME.message);
end

%% ---------- Signaux collés & erreurs (tous les intermédiaires) ---------------
try
    needFEM = (exist('EITFEM','class') == 8) && ~isempty(S) && isfield(S,'E') && isfield(S,'params');
    hasData = isfield(R,'Umeas') && ~isempty(R.Umeas) && isfield(R,'Ufwd') && ~isempty(R.Ufwd);
    if needFEM && hasData
        Ne = S.params.Ne;
        K  = floor(numel(R.Umeas)/Ne);

        fem = EITFEM(ForwardMesh1st(g*1e-3, H, S.E)); fem.mode='current';
        fem.zeta = S.params.z_contact*ones(Ne,1);

        if isfield(S,'Ipat') && ~isempty(S.Ipat), base=S.Ipat;
        else
            if exist('EITSim','class') == 8 || exist('EITSim','class') == 2
                base=EITSim.buildTrigPattern(Ne, S.params.I_amp);
            else
                base=inj_pattern_trig(Ne, S.params.I_amp);
            end
        end
        if size(base,2) < K, base=repmat(base,1,ceil(K/size(base,2))); end
        inj = base(:,1:K); fem.Iel = inj(:);

        U_meas = R.Umeas(:);
        U_rec  = R.Ufwd(:);
        scaleU = 1; if isfield(R,'scaleU') && ~isempty(R.scaleU), scaleU = R.scaleU; end

        fem.Uel = U_meas; fem.SetInvGamma(5e-6,8e-3);
        sigma0_nodal = [];
        if isfield(R,'sigma0') && ~isempty(R.sigma0), sigma0_nodal = R.sigma0(:); end
        if ~isempty(sigma0_nodal)
            U_init  = fwd_norm_from_sigma(fem, sigma0_nodal, Ne, scaleU);
        else
            U_init = R.Ufwd(:);
        end

        % Intermédiaires (tous)
        U_inters = []; labIters = [];
        if exist(iterDir,'dir')
            D = dir(fullfile(iterDir,'iter_*.mat'));
            [~,ord] = sort({D.name}); D = D(ord);
            U_inters = zeros(Ne*K, numel(D));
            labIters = zeros(numel(D),1);
            for i=1:numel(D)
                Si = load(fullfile(iterDir, D(i).name));
                sig_i = [];
                if isfield(Si,'sigma_iter_nodal') && ~isempty(Si.sigma_iter_nodal)
                    sig_i = Si.sigma_iter_nodal(:);
                elseif isfield(Si,'sigma_iter_tri') && ~isempty(Si.sigma_iter_tri)
                    sig_i = nodal_from_tri(g,H,Si.sigma_iter_tri(:));
                end
                if isempty(sig_i), continue; end
                U_inters(:,i) = fwd_norm_from_sigma(fem, sig_i, Ne, scaleU);
                if isfield(Si,'iter_global') && ~isempty(Si.iter_global)
                    labIters(i) = Si.iter_global;
                else
                    labIters(i) = i;
                end
            end
            keep = any(U_inters~=0,1);
            U_inters = U_inters(:,keep);
            labIters = labIters(keep);
        end

        x = 1:(Ne*K);
        fGlue = figure('Color','w','Name','Signals','NumberTitle','off');
        plot(x, U_meas, '-',  'DisplayName','measured'); hold on;
        plot(x, U_init, '-',  'DisplayName','init');
        plot(x, U_rec,  '-',  'DisplayName','reconstruction');
        if ~isempty(U_inters)
            for i=1:size(U_inters,2)
                plot(x, U_inters(:,i), '-', 'HandleVisibility','off', 'Color',[0.5 0.5 0.5 0.35]);
            end
        end
        for m = 1:K-1, xline(m*Ne, ':', 'HandleVisibility','off'); end
        grid on; xlim([1 Ne*K]);
        xlabel('Electrode index across all patterns (1..Ne*K)'); ylabel('V (normalized)');
        title(sprintf('Signals (Ne=%d, K=%d) — %s slice %d', Ne, K, patient_id, z_slice));
        legend('Location','best');
        if SAVE_PNG
            exportgraphics(fGlue, fullfile(plotDir, ...
                sprintf('signals_slice%03d_initConst_%s.png', z_slice, num2tag(INIT_CONST))), ...
                'Resolution', DPI);
        end

        if ~isempty(U_inters)
            fGlueOnly = figure('Color','w','Name','Signals (intermediates only)','NumberTitle','off');
            hold on;
            nI = size(U_inters,2);
            C = hsv(max(nI,2));
            for i=1:nI
                plot(x, U_inters(:,i), '-', 'Color', C(i,:), 'DisplayName', sprintf('iter %d', labIters(i)));
            end
            for m = 1:K-1, xline(m*Ne, ':', 'HandleVisibility','off'); end
            grid on; xlim([1 Ne*K]);
            xlabel('Electrode index across all patterns (1..Ne*K)'); ylabel('V (normalized)');
            title(sprintf('Signals — intermediates only (n=%d)', nI));
            if nI <= 15, legend('Location','eastoutside'); end
            if SAVE_PNG
                exportgraphics(fGlueOnly, fullfile(plotDir, ...
                    sprintf('signals_intermediates_only_slice%03d_initConst_%s.png', z_slice, num2tag(INIT_CONST))), ...
                    'Resolution', DPI);
            end
        end

        err_init = abs(U_init - U_meas);
        err_rec  = abs(U_rec  - U_meas);
        fErr = figure('Color','w','Name','Signal errors','NumberTitle','off');
        plot(x, err_init, '-', 'DisplayName','|init - meas|'); hold on;
        plot(x, err_rec,  '-',  'DisplayName','|reco - meas|');
        if ~isempty(U_inters)
            for i=1:size(U_inters,2)
                ei = abs(U_inters(:,i) - U_meas);
                plot(x, ei, '-', 'HandleVisibility','off', 'Color',[0.6 0.6 0.6 0.25]);
            end
        end
        for m = 1:K-1, xline(m*Ne, ':', 'HandleVisibility','off'); end
        grid on; xlim([1 Ne*K]);
        xlabel('Electrode index across all patterns (1..Ne*K)'); ylabel('abs error (normalized)');
        ttl = sprintf('Signal errors (RMSE init=%.3g, reco=%.3g)', ...
                      sqrt(mean(err_init.^2)), sqrt(mean(err_rec.^2)));
        title(ttl);
        legend('Location','best');
        if SAVE_PNG
            exportgraphics(fErr, fullfile(plotDir, ...
                sprintf('signal_errors_slice%03d_initConst_%s.png', z_slice, num2tag(INIT_CONST))), ...
                'Resolution', DPI);
        end
    end
catch ME
    warning('Plots error non générés: %s', ME.message);
end

%% ================= helpers locaux =================
function tri_colorplot(g,H,tri_vals,clip_range,displayMode)
    if size(g,2)~=2, error('tri_colorplot: g doit être [Nnodes x 2] (mm).'); end
    if size(H,2)~=3, error('tri_colorplot: H doit être [Ntri x 3].'); end
    tri_vals = tri_vals(:);
    if numel(tri_vals)~=size(H,1)
        error('tri_colorplot: tri_vals doit avoir Ntri (=size(H,1)) éléments.');
    end
    patch('Faces',H,'Vertices',g, ...
          'FaceVertexCData', tri_vals, ...
          'FaceColor','flat', ...
          'EdgeColor','none');
    axis equal tight;
    if exist('displayMode','var') && ischar(displayMode) && strcmpi(displayMode,'radiological')
        set(gca,'XDir','reverse');
    else
        set(gca,'XDir','normal');
    end
    if exist('turbo','file'), colormap(turbo); else, colormap(parula); end
    if ~isempty(clip_range) && numel(clip_range)==2 && all(isfinite(clip_range)), caxis(clip_range); end
    set(gca,'XTick',[],'YTick',[]);
    xlabel('x (mm)'); ylabel('y (mm)');
end

function tri = tri_from_nodal_mode(H, sigma_nod, levels)
    if nargin<3 || isempty(levels)
        lv = unique(sigma_nod(:));
        if numel(lv)>20, levels = linspace(min(sigma_nod), max(sigma_nod), 5);
        else,             levels = lv(:).'; end
    end
    L = levels(:).';
    [~, lab_nod] = min(abs(sigma_nod(:) - L), [], 2);
    tri_lab = mode([lab_nod(H(:,1)) lab_nod(H(:,2)) lab_nod(H(:,3))], 2);
    tri = L(tri_lab);
end

function levels = infer_levels_from_packs(R, S)
    levels = [];
    if ~isempty(S) && isfield(S,'params') && isfield(S.params,'cond')
        c = S.params.cond; fn = fieldnames(c); vals = zeros(numel(fn),1);
        for i=1:numel(fn), vals(i)=c.(fn{i}); end
        levels = unique(sort(vals));
    elseif isfield(R,'initInfo') && isfield(R.initInfo,'levels') && ~isempty(R.initInfo.levels)
        levels = R.initInfo.levels(:);
    end
    if isempty(levels)
        s = R.sigma0(:);
        qs = quantile(s, [0.2 0.4 0.6 0.8]);
        levels = unique(sort([min(s); qs(:); max(s)]));
    end
end

function U_norm = fwd_norm_from_sigma(fem, sigma_nodal_vec, Ne, scaleU)
    Uraw = fem.SolveForwardVec(sigma_nodal_vec(:));
    K = numel(Uraw)/Ne;
    V = reshape(Uraw, Ne, K);
    V = V - mean(V,1);
    scale = max(scaleU, 1);
    U_norm = (V(:)) ./ scale;
end

function nodal = nodal_from_tri(g,H,triVals)
    N = size(g,1);
    nodal = zeros(N,1);
    cnt   = zeros(N,1);
    triVals = triVals(:);
    for t=1:size(H,1)
        v = triVals(t);
        idx = H(t,:).';
        nodal(idx) = nodal(idx) + v;
        cnt(idx)   = cnt(idx) + 1;
    end
    cnt(cnt==0) = 1;
    nodal = nodal ./ cnt;
end

function inj = inj_pattern_trig(Ne, I_amp)
    inj = zeros(Ne, Ne);
    idx = (0:Ne-1)';
    for k = 1:Ne
        phase = 2*pi*(k-1)/Ne;
        v = sin( 2*pi*idx/Ne + phase );
        v = v - mean(v);
        if max(abs(v)) > 0, v = v / max(abs(v)); end
        v = v - mean(v);
        inj(:,k) = I_amp * v;
    end
end
