% ============================
% File: Problem_Forward.m
% ============================
% PROBLÈME DIRECT EIT (batch-only)
% -------------------------------------------------------------------------
% BUT
%   Préparer et exécuter une simulation FORWARD EIT pour une coupe CT :
%     1) Extraire les contours/organes en coordonnées mm
%     2) Construire un maillage 2D (PDE Toolbox)
%     3) Poser des électrodes régulièrement sur le bord
%     4) Définir des paramètres EIT (Ne, I_amp, z_contact, conductivités)
%     5) Simuler le forward (OOEIT) et sauvegarder les .mat nécessaires
%
% IMPORTANT
%   - AUCUNE FIGURE n'est affichée ici (batch). Les visualisations se font
%     dans Plot_Problem_Forward.m.
%   - Les chemins d'entrée/sortie par défaut ciblent Outputs/<patient>/slice_xxx.
%
% PRÉREQUIS / DÉPENDANCES
%   - Dossier 'src' contenant les utilitaires (extract_slice_polys, MeshBuilder,
%     buildElectrodesFromContour, EITSim, etc.)
%   - OOEIT (EITFEM, ForwardMesh1st, …) accessible dans le path MATLAB.
%   - Données NIfTI : Data_set/<patient_id>/ct.nii.gz + segmentations/*.nii(.gz)
%
% SORTIES (dans Outputs/<patient_id>/slice_xxx/)
%   - mesh/mesh_sliceXXX.mat        : g, H, triGroup, shapes, domain, contour
%   - E_electrodes.mat              : E (cell d'arêtes), el_centers
%   - cond_values.mat               : cond (struct des conductivités par groupe)
%   - viz_config.mat                : limites et options de rendu
%   - Imeas.mat, U_cell.mat, Vmat.mat, Ipat.mat, Efield_tri.mat
%   - eit_pack.mat                  : pack complet standard (via EITSim.save_results)
%   - data_sources.mat              : méta (ct_file, seg_dir, z_slice, patient_id)
%
% PIÈGES COURANTS
%   - Vérifier que 'src' et OOEIT sont dans le path (addpath ci-dessous).
%   - displayMode = 'radiological' inverse l'axe X (convention DICOM).
%   - Adapter 'mesh_params.targetSize' et la largeur d'électrode 'el_width'
%     selon l'échelle (mm) et le périmètre du contour.
% -------------------------------------------------------------------------

clear; close all; clc;

% --- 0) PATHS & MODE SILENCIEUX ------------------------------------------
addpath('src', genpath('src'));  % utilitaires du dépôt (FFT lissage, NIfTI, mesh, etc.)
addpath(genpath('/Users/anis/Documents/StageInria/Code/OOEIT-main')); % OOEIT

% Aucune figure ne doit s'afficher dans ce script
set(groot,'defaultFigureVisible','off');

% --- 1) PARAMÈTRES UTILISATEUR -------------------------------------------
patient_id  = 's0011';           % ID du patient (doit exister dans Data_set/)
z_slice     = 301;               % indice de coupe
displayMode = 'neurological';    % ou 'radiological' (X inversé)

% Dossiers de sortie (créés si absents)
rootOut     = fullfile('Outputs', patient_id, sprintf('slice_%03d', z_slice));
meshDir     = fullfile(rootOut, 'mesh');
if ~exist(rootOut,'dir'),  mkdir(rootOut);  end
if ~exist(meshDir,'dir'),  mkdir(meshDir);  end

% --- 2) EXTRACTION CONTOURS / SHAPES -------------------------------------
% Produit :
%   - domain     : polyshape du "soft tissue" (extérieur - organes)
%   - shapes     : struct de polyshape par groupes (Heart, Lung, …)
%   - ext_smooth : contour extérieur lissé [N x 2] en mm (fermé)
%   - info       : metadata NIfTI du CT
[domain, shapes, ext_smooth, info] = extract_slice_polys(patient_id, z_slice); %#ok<ASGLU>

% --- 3) CONSTRUCTION DU MAILLAGE -----------------------------------------
% Paramètres de maillage (PDE Toolbox)
mesh_params = struct('targetSize', 5, ... % taille max des triangles (mm)
                     'minArea_mm2', 300, ...
                     'Hgrad', 1.4);      % lissage de la taille locale

% Bâtit la géométrie PDE (contour ext + organes), triangule, étiquette par groupe
meshObj = MeshBuilder.from_polys(shapes, domain, ext_smooth, z_slice, mesh_params);

% Sauvegarde du mesh pour les scripts de plotting/inverse
% -> écrit mesh/mesh_sliceXXX.mat avec g, H, triGroup, shapes, domain, contour
meshObj.saveMesh(meshDir, z_slice);

% --- 4) ÉLECTRODES --------------------------------------------------------
% Place Ne électrodes d'arc el_width (mm) sur le bord libre du maillage
Ne = 16;
el_width_mm = 33;
[E, el_centers] = buildElectrodesFromContour(meshObj.g, meshObj.contour, Ne, el_width_mm, meshObj.H);
save(fullfile(rootOut, 'E_electrodes.mat'), 'E', 'el_centers');

% --- 5) PARAMÈTRES EIT ----------------------------------------------------
% Conductivités par groupes (S/m) + autres paramètres forward
eit_params = struct( ...
    'Ne',        Ne, ...
    'el_width',  el_width_mm, ...     % mm
    'I_amp',     0.35e-3, ...         % A
    'z_contact', 0.05, ...            % Ohm·m^2
    'groups',    {meshObj.groups}, ...
    'cond',      struct( ...
        'SoftTissue',            0.35, ... % ~0.25–0.50 (glandes, moelle, etc.)
        'Heart',                 0.50, ...
        'Lung',                  0.15, ...
        'Trachea',               0.12, ... % approx. cartilage
        'Bone',                  0.03, ... % ~0.02–0.05
        'Cartilage',             0.12, ... % ~0.08–0.15
        'Blood',                 0.70, ...
        'LiverSpleenPancreas',   0.18, ... % ~0.15–0.20
        'Kidney',                0.35, ...
        'Other',                 0.30 ...
    ), ...
    'noiseRel',  1e-3, ...            % bruit relatif pour la synthèse
    'rngSeed',   123 ...              % seed pour reproductibilité
);

% Sauvegarde pratique des conductivités (utilisé côté inverse)
cond = eit_params.cond; 
save(fullfile(rootOut, 'cond_values.mat'), 'cond');

% --- 6) BORNES / VIZ CONFIG GLOBALES -------------------------------------
% Limites de colorbar basées sur les valeurs cond (pour plotting)
vals = struct2array(eit_params.cond);
SIGMA_BOUNDS = [min(vals), max(vals)];

viz_config = struct('displayMode', displayMode, ...
                    'sigma_clim',  SIGMA_BOUNDS, ...
                    'cmap',        'turbo');
save(fullfile(rootOut,'viz_config.mat'),'viz_config');

% --- 7) SIMULATION FORWARD (AUCUN PLOT) ----------------------------------
% Instancie le simulateur (convertit mm->m pour OOEIT, place Iel, etc.)
eitSim = EITSim(meshObj.g, meshObj.H, meshObj.triGroup, ...
                meshObj.domain, eit_params, E, el_centers);

% Calcule :
%   sigma_tri, sigma_used, Imeas_clean, Imeas, U_cell (par motif),
%   Vmat (Ne x K), Ipat (Ne x K), E_mag_tri, gradUx_tri, gradUy_tri
eitSim = eitSim.simulate_forward_full();

% Sauvegardes unitaires (faciles à charger)
Imeas = eitSim.Imeas;                        %#ok<NASGU>
U_cell = eitSim.U_cell;                      %#ok<NASGU>
Vmat   = eitSim.Vmat;                        %#ok<NASGU>
Ipat   = eitSim.Ipat;                        %#ok<NASGU>
E_mag_tri  = eitSim.E_mag_tri;               %#ok<NASGU>
gradUx_tri = eitSim.gradUx_tri;              %#ok<NASGU>
gradUy_tri = eitSim.gradUy_tri;              %#ok<NASGU>
save(fullfile(rootOut, 'Imeas.mat'),      'Imeas');
save(fullfile(rootOut, 'U_cell.mat'),     'U_cell');
save(fullfile(rootOut, 'Vmat.mat'),       'Vmat');
save(fullfile(rootOut, 'Ipat.mat'),       'Ipat');
save(fullfile(rootOut, 'Efield_tri.mat'), 'E_mag_tri','gradUx_tri','gradUy_tri');

% Pack standard complet (pour Plot/Inverse)
eitSim.save_results(rootOut);

% --- 8) CONTEXTE CT/SEGMENTATIONS ----------------------------------------
% Chemins vers les sources brutes (utile pour traçabilité et scripts plotting)
ct_file = fullfile('Data_set', patient_id, 'ct.nii.gz');
seg_dir = fullfile('Data_set', patient_id, 'segmentations');
meta = struct('ct_file', ct_file, 'seg_dir', seg_dir, 'z_slice', z_slice, ...
              'patient_id', patient_id);
save(fullfile(rootOut, 'data_sources.mat'), 'meta');

% Restaure l'affichage par défaut pour la session suivante
set(groot,'defaultFigureVisible','on');

disp('Problem_Forward terminé (batch-only, aucun plot affiché).');
