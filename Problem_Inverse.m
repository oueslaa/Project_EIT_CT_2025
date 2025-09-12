% ============================
% File: Problem_Inverse.m
% ============================
% RECONSTRUCTION EIT (inverse problem) — slice unique
%   - Transplante les organes de slice_prev dans le contour ext. de slice_this
%   - Construit un maillage INVERSE propre (contour seul)
%   - Crée une init σ0 structurée (discrète) depuis les polygones (pas constante)
%   - Lance NOSER/GN sur le maillage inverse (découplé du forward)

clear; close all; clc;

% --- PATHS ---------------------------------------------------------------
addpath(genpath('src_anis'));
addpath(genpath('/Users/anis/Documents/StageInria/Code/OOEIT-main')); % <- adapte

set(groot,'defaultFigureVisible','off');  % pas de figures ici

%% -------- Cas -----------------------------------------------------------
patient_this = 's0011';  slice_this = 301;   % cible
patient_prev = 's0011';  slice_prev = 310;   % donneur d'organes (contours)

%% -------- Hyperparamètres « best » --------------------------------------
BEST_MU    = 0.05;    % NOSER mu
BEST_ALPHA = 1e-3;    % alpha L2

N_NOSER = 1;          % nb étapes NOSER
N_GN    = 3;          % nb itérations GN

%% -------- Dossiers d'E/S ------------------------------------------------
rootThis = fullfile('Outputs',patient_this,sprintf('slice_%03d',slice_this));
rootPrev = fullfile('Outputs',patient_prev,sprintf('slice_%03d',slice_prev));

if ~exist(rootThis,'dir'), mkdir(rootThis); end

packFile = fullfile(rootThis,'eit_pack.mat');    % pack FORWARD existant (mesh riche)
meshFwd  = fullfile(rootThis,'mesh',sprintf('mesh_slice%03d.mat',slice_this)); % mesh riche cible (contient contour)
meshPrev = fullfile(rootPrev,'mesh',sprintf('mesh_slice%03d.mat',slice_prev)); % mesh donneur (contient shapes)

assert(isfile(packFile), 'Pack forward cible manquant: %s', packFile);
assert(isfile(meshFwd),  'Mesh forward cible manquant: %s', meshFwd);
assert(isfile(meshPrev), 'Mesh donneur manquant: %s', meshPrev);

% Répertoires inverse
invRoot = fullfile(rootThis,'inverse',sprintf('prev_%s_%03d',patient_prev,slice_prev));
plotDir = fullfile(invRoot,'plots_inverse');
iterDir = fullfile(plotDir,'iterations');
if ~exist(invRoot,'dir'), mkdir(invRoot); end
if ~exist(plotDir,'dir'), mkdir(plotDir); end
if ~exist(iterDir,'dir'), mkdir(iterDir); end

recFile = fullfile(invRoot,'reconstruction_inverse.mat');

%% -------- Clip dynamique depuis pack -----------------------------------
S = load(packFile);
vals = struct2array(S.params.cond);
lo = max(1e-3, min(vals)); hi = max(vals); span = hi - lo;
CLIP = [max(1e-3, lo-0.02*span), hi+0.02*span];

%% -------- Options solveur -----------------------------------------------
opts = struct();

% --- Re-maillage inverse (contour seul) + transplantation d'organes prev ---
opts.inverse_from_prev_organs = true;     % indicateur (ok si non lu par simulate)
opts.inv_Hmax  = 5;                       % mm
opts.inv_Hgrad = 1.3;
opts.rebuild_inverse_mesh = false;

% --- Warp contours prev->this (TPS) ---
opts.doWarp = true;
opts.tps    = struct('enable',true,'lambda',1e-3);

% --- Init MORPHOLOGIQUE (pas constante !) ---
opts.init = struct( ...
    'constant_enable', false, ...   % <<< IMPORTANT: désactive init homogène
    'discrete',        true, ...
    'knn_k',           11, ...
    'clean_iters',     3, ...
    'band_mm',         3.0, ...
    'island_min_area_mm2', 200, ...
    'smooth_iters',    0, ...
    'smooth_lambda',   0, ...
    'calib_enable',    false, ...
    'calib_clip',      [1 1], ...
    'noser_warmstart', struct('enable',false) );

% --- Paramétrisation / régularisation ---
opts.param_by   = 'node';                   % 'node' ou 'tri'
opts.reg_alpha  = (N_GN>0) * BEST_ALPHA;    % L2 sur θ si GN>0
opts.reg_on     = 'theta';

% --- Données / bruit ---
opts.clip       = CLIP;
opts.constNoise = 3e-5;
opts.relNoise   = 1.2e-2;

% --- Itérations ---
opts.iters_per_stage = max(0, N_GN);
opts.noser_only = struct( ...
    'enable',       N_NOSER>0, ...
    'iters',        max(0, N_NOSER), ...
    'lambda',       BEST_MU, ...
    'lambda_mode',  'scaled-median', ...
    'use_diag',     true, ...
    'max_rel_step', 0.20, ...
    'linesearch',   false, ...
    'alpha_l2',     BEST_ALPHA );

% --- Snapshots ---
opts.snapshots_enable = true;
opts.cleanIterDir     = true;

fprintf('[INVERSE] target=(%s,%03d) | prev=(%s,%03d)\n', patient_this, slice_this, patient_prev, slice_prev);
fprintf('[INVERSE] Using tuned best: mu=%.6g, alpha=%.6g | clip=[%.4g, %.4g]\n', ...
        BEST_MU, BEST_ALPHA, CLIP(1), CLIP(2));
fprintf('[INVERSE] Plan: NOSER=%d step(s), GN=%d it(s)\n', N_NOSER, N_GN);
fprintf('[INVERSE] Output dir: %s\n', invRoot);

%% -------- Lancement -----------------------------------------------------
rec = simulate_inverse_slice(patient_this, slice_this, patient_prev, slice_prev, opts, iterDir);

%% -------- Sauvegarde ----------------------------------------------------
save(recFile, '-struct', 'rec');
mis = NaN; if isfield(rec,'misfit') && ~isempty(rec.misfit), mis = rec.misfit; end
fprintf('\n[INVERSE] Reconstruction done. Relative misfit = %.6g\nFile: %s\n', mis, recFile);
