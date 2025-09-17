function Problem_Forward_fn(patient_id, z_slice, varargin)
%PROBLEM_FORWARD_FN  Exécute le forward EIT pour (patient_id, z_slice).
% Usage:
%   Problem_Forward_fn('s0011', 320)
%
% Options (Name-Value) :
%   'DisplayMode'  : 'neurological' (def) ou 'radiological'
%   'Ne'           : nombre d'électrodes (def 16)
%   'ElWidthMM'    : largeur (mm) d'une électrode (def 33)
%   'Iamp'         : amplitude courant (A) (def 0.35e-3)
%   'Zcontact'     : Ohm·m^2 (def 0.05)
%   'NoiseRel'     : bruit relatif (def 1e-3)
%   'RngSeed'      : entier ou [] (def 123)

p = inputParser;
addParameter(p,'DisplayMode','neurological');
addParameter(p,'Ne',16);
addParameter(p,'ElWidthMM',33);
addParameter(p,'Iamp',0.35e-3);
addParameter(p,'Zcontact',0.05);
addParameter(p,'NoiseRel',1e-3);
addParameter(p,'RngSeed',123);
parse(p,varargin{:});
opt = p.Results;

% --- PATHS & MODE SILENCIEUX ---------------------------------------------
addpath('src_anis', genpath('src_anis'));   % utilitaires
ooeitDir = fullfile(fileparts(pwd), 'OOEIT-main');
if exist(ooeitDir,'dir'), addpath(genpath(ooeitDir)); end

set(groot,'defaultFigureVisible','off');

% --- Sorties --------------------------------------------------------------
rootOut = fullfile('Outputs', patient_id, sprintf('slice_%03d', z_slice));
meshDir = fullfile(rootOut, 'mesh');
if ~exist(rootOut, 'dir'),  mkdir(rootOut);  end
if ~exist(meshDir, 'dir'),  mkdir(meshDir);  end

% --- 2) Extraction contours / shapes -------------------------------------
[domain, shapes, ext_smooth, ~] = extract_slice_polys(patient_id, z_slice); %#ok<ASGLU>

% --- 3) Maillage ----------------------------------------------------------
mesh_params = struct('targetSize',5, 'minArea_mm2',300, 'Hgrad',1.4);
meshObj = MeshBuilder.from_polys(shapes, domain, ext_smooth, z_slice, mesh_params);
meshObj.saveMesh(meshDir, z_slice);

% --- 4) Électrodes --------------------------------------------------------
Ne = opt.Ne; el_width_mm = opt.ElWidthMM;
[E, el_centers] = buildElectrodesFromContour(meshObj.g, meshObj.contour, Ne, el_width_mm, meshObj.H);
save(fullfile(rootOut,'E_electrodes.mat'), 'E','el_centers');

% --- 5) Paramètres EIT ----------------------------------------------------
eit_params = struct( ...
    'Ne',        Ne, ...
    'el_width',  el_width_mm, ...
    'I_amp',     opt.Iamp, ...
    'z_contact', opt.Zcontact, ...
    'groups',    {meshObj.groups}, ...
    'cond',      struct( ...
        'SoftTissue',            0.35, ...
        'Heart',                 0.50, ...
        'Lung',                  0.15, ...
        'Trachea',               0.12, ...
        'Bone',                  0.03, ...
        'Cartilage',             0.12, ...
        'Blood',                 0.70, ...
        'LiverSpleenPancreas',   0.18, ...
        'Kidney',                0.35, ...
        'Other',                 0.30 ...
    ), ...
    'noiseRel',  opt.NoiseRel, ...
    'rngSeed',   opt.RngSeed ...
);

cond = eit_params.cond; 
save(fullfile(rootOut,'cond_values.mat'),'cond');

% --- 6) Viz config --------------------------------------------------------
vals = struct2array(eit_params.cond);
viz_config = struct('displayMode', opt.DisplayMode, ...
                    'sigma_clim',  [min(vals), max(vals)], ...
                    'cmap',        'turbo');
save(fullfile(rootOut,'viz_config.mat'),'viz_config');

% --- 7) Simulation forward ------------------------------------------------
eitSim = EITSim(meshObj.g, meshObj.H, meshObj.triGroup, meshObj.domain, ...
                eit_params, E, el_centers);
eitSim = eitSim.simulate_forward_full();

Imeas = eitSim.Imeas;            %#ok<NASGU>
U_cell = eitSim.U_cell;          %#ok<NASGU>
Vmat   = eitSim.Vmat;            %#ok<NASGU>
Ipat   = eitSim.Ipat;            %#ok<NASGU>
E_mag_tri  = eitSim.E_mag_tri;   %#ok<NASGU>
gradUx_tri = eitSim.gradUx_tri;  %#ok<NASGU>
gradUy_tri = eitSim.gradUy_tri;  %#ok<NASGU>
save(fullfile(rootOut,'Imeas.mat'),      'Imeas');
save(fullfile(rootOut,'U_cell.mat'),     'U_cell');
save(fullfile(rootOut,'Vmat.mat'),       'Vmat');
save(fullfile(rootOut,'Ipat.mat'),       'Ipat');
save(fullfile(rootOut,'Efield_tri.mat'), 'E_mag_tri','gradUx_tri','gradUy_tri');

eitSim.save_results(rootOut);

% --- 8) Sources -----------------------------------------------------------
ct_file = fullfile('Data_set', patient_id, 'ct.nii.gz');
if ~isfile(ct_file), ct_file = fullfile('Data_set', patient_id, 'ct.nii'); end
seg_dir = fullfile('Data_set', patient_id, 'segmentations');
meta = struct('ct_file', ct_file, 'seg_dir', seg_dir, 'z_slice', z_slice, 'patient_id', patient_id);
save(fullfile(rootOut,'data_sources.mat'),'meta');

set(groot,'defaultFigureVisible','on');
end
