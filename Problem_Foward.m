% ============================
% File: Problem_Foward.m
% ============================
% Problème direct EIT : prépare toutes les données et enregistre
% uniquement des .mat nécessaires. AUCUNE FENÊTRE N'EST AFFICHÉE.
% Les visualisations sont faites dans Plot_Problem_Foward.m

clear; close all; clc;
addpath('src', genpath('src'));
addpath(genpath('/Users/anis/Documents/StageInria/Code/OOEIT-main'));

% ---- Important : aucune figure ne doit s'afficher ici
set(groot,'defaultFigureVisible','off');

% ---------- Paramètres utilisateur ----------
patient_id  = 's0011';
z_slice     = 202;
displayMode = 'neurological';   % ou 'radiological'

rootOut     = fullfile('Outputs', patient_id, sprintf('slice_%03d', z_slice));
meshDir     = fullfile(rootOut, 'mesh');
if ~exist(rootOut,'dir'), mkdir(rootOut); end
if ~exist(meshDir,'dir'), mkdir(meshDir); end

% ---------- 1) Extraction contours / shapes  ----------
[domain, shapes, ext_smooth, info] = extract_slice_polys(patient_id, z_slice); %#ok<ASGLU>

% ---------- 2) Construction du maillage  ----------
mesh_params = struct('targetSize', 5, 'minArea_mm2', 300, 'Hgrad', 1.4);
meshObj = MeshBuilder.from_polys(shapes, domain, ext_smooth, z_slice, mesh_params);

% Sauvegarde du mesh pour réutilisation par les scripts Plot/Inverse
meshObj.saveMesh(meshDir, z_slice);
% (Écrit : mesh/mesh_sliceXXX.mat avec g, H, triGroup, shapes, domain, contour)

% ---------- 3) Électrodes  ----------
[E, el_centers] = buildElectrodesFromContour(meshObj.g, meshObj.contour, 16, 33, meshObj.H);
save(fullfile(rootOut, 'E_electrodes.mat'), 'E', 'el_centers');

% ---------- 4) Paramètres EIT (réutilisables) ----------
eit_params = struct( ...
    'Ne',        16, ...
    'el_width',  33, ...         % mm
    'I_amp',     0.35e-3, ...    % A
    'z_contact', 0.05, ...       % Ohm·m^2
    'groups',    {meshObj.groups}, ...
    'cond',      struct('SoftTissue',0.30,'Heart',0.50,'Lung',0.15,'Trachea',0.15,'Bone',0.05,'Other',0.28), ...
    'noiseRel',  1e-3, ...
    'rngSeed',   123 ...
);

% Sauvegarde des conductivités par organe (utile pour l'inverse)
cond = eit_params.cond; 
save(fullfile(rootOut, 'cond_values.mat'), 'cond');

% ---------- 5) Bornes/colormap à réutiliser pour TOUS les plots ----------
SIGMA_BOUNDS = [ ...
    min([eit_params.cond.SoftTissue, eit_params.cond.Bone, eit_params.cond.Lung, eit_params.cond.Trachea, eit_params.cond.Other]), ...
    max([eit_params.cond.SoftTissue, eit_params.cond.Heart, eit_params.cond.Other]) ...
];
viz_config = struct('displayMode', displayMode, 'sigma_clim', SIGMA_BOUNDS, 'cmap', 'turbo');
save(fullfile(rootOut,'viz_config.mat'),'viz_config');

% ---------- 6) Simulation FORWARD (pas de plots) ----------
eitSim = EITSim(meshObj.g, meshObj.H, meshObj.triGroup, meshObj.domain, eit_params, E, el_centers);
eitSim = eitSim.simulate_forward();

% Sauvegardes (aucun PNG créé ici)
Imeas = eitSim.Imeas; 
save(fullfile(rootOut, 'Imeas.mat'), 'Imeas');

% Pack standard complet pour réutilisation (Plot/Inverse)
% (g,H,triGroup,E,el_centers,sigma_tri,params,Imeas,domain)
eitSim.save_results(rootOut);

% ---------- 7) Contexte CT/segmentation (pour le script Plot) ----------
ct_file = fullfile('Data_set', patient_id, 'ct.nii.gz');
seg_dir = fullfile('Data_set', patient_id, 'segmentations');
meta = struct('ct_file', ct_file, 'seg_dir', seg_dir, 'z_slice', z_slice, ...
              'patient_id', patient_id);
save(fullfile(rootOut, 'data_sources.mat'), 'meta');

disp('Problem_Foward terminé (batch-only, aucun plot affiché).');
