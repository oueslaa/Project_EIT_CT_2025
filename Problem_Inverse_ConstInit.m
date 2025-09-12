% ======================================
% File: Problem_Inverse_ConstInit.m
% ======================================
% Reconstruction avec init homogène (sigma constant) + hyperparams OPTIMAUX.
% Sauvegarde:
%   Outputs/<patient>/slice_<ZZZ>/inverse/initConst_<val>/

clear; close all; clc;
addpath(genpath('src_anis'));
addpath(genpath('/Users/anis/Documents/StageInria/Code/OOEIT-main'));  % adapte si besoin
set(groot,'defaultFigureVisible','off');  % aucun plot ici

%% -------- Choix cas --------
patient_this = 's0011';
slice_this   = 301;

% Pour init constante on peut ignorer prev, mais la signature l'exige :
patient_prev = patient_this;
slice_prev   = slice_this;

%% -------- Hyperparamètres OPTIMAUX (depuis tuning) --------
BEST_MU    = 0.05;   % μ* pour NOSER
BEST_ALPHA = 1e-3;   % α* (ridge NOSER + régularisation GN)

% Planning d'itérations (si tu veux pousser après init)
N_NOSER = 1;
N_GN    = 3;

%% -------- Chemins & pack --------
rootThis = fullfile('Outputs',patient_this,sprintf('slice_%03d',slice_this));  % pack/mesh
packFile = fullfile(rootThis,'eit_pack.mat');
assert(isfile(packFile), 'Pack forward introuvable: %s', packFile);

S = load(packFile);

% Clip dynamique robuste depuis conds du pack
vals = struct2array(S.params.cond);
lo = max(1e-3, min(vals)); hi = max(vals); span = hi - lo;
CLIP = [max(1e-3, lo - 0.02*span), hi + 0.02*span];

% Valeur pour l'init constante
if isfield(S,'params') && isfield(S.params,'cond') && isfield(S.params.cond,'SoftTissue')
    INIT_CONST = S.params.cond.SoftTissue;
elseif ~isempty(vals)
    INIT_CONST = median(vals);
else
    INIT_CONST = mean(CLIP);
end

% Dossiers de sortie
num2tag = @(x) strrep(sprintf('%.4g',x),'.','p');
invRoot = fullfile(rootThis,'inverse', sprintf('initConst_%s', num2tag(INIT_CONST)));
plotDir = fullfile(invRoot,'plots_inverse');
iterDir = fullfile(plotDir,'iterations');
if ~exist(invRoot,'dir'), mkdir(invRoot); end
if ~exist(plotDir,'dir'), mkdir(plotDir); end
% purge snapshots avant run
if exist(iterDir,'dir'), try, rmdir(iterDir,'s'); catch, end, end
mkdir(iterDir);

recFile = fullfile(invRoot,'reconstruction_inverse.mat');

fprintf('[INVERSE-CONST] target=(%s,%03d)\n', patient_this, slice_this);
fprintf('[INVERSE-CONST] init_const=%.6g | mu=%.6g | alpha=%.6g | clip=[%.4g, %.4g]\n', ...
        INIT_CONST, BEST_MU, BEST_ALPHA, CLIP(1), CLIP(2));
fprintf('[INVERSE-CONST] plan: NOSER=%d, GN=%d | outdir=%s\n', N_NOSER, N_GN, invRoot);

%% ====== Options pour simulate_inverse_slice =============================
opts = struct();

% --- IMPORTANT : pas d'organes, maillage inverse neutre (contour seul) ---
opts.inverse_from_prev_organs = false;     % <<< désactive la transplantation d'organes
opts.inv_Hmax  = 5;
opts.inv_Hgrad = 1.3;
opts.rebuild_inverse_mesh = false;

% --- Paramétrisation / régularisation ---
opts.param_by   = 'node';                   % 'node' ou 'tri'
opts.reg_alpha  = (N_GN>0) * BEST_ALPHA;

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

% --- INIT HOMOGÈNE (fait foi) ---
opts.init = struct( ...
    'constant_enable', true, ...      % <- utilisé dans simulate_inverse_slice
    'soft_value',      INIT_CONST, ...% <- valeur homogène
    'discrete',        true, ...
    'knn_k',           0, ...
    'quant',           struct('clean_iters',0,'centers',[]), ...
    'auto_blend',      false, 'blend_max',0, 'blend_homog',0, ...
    'smooth_iters',    0, 'smooth_lambda',0, ...
    'calib_enable',    false, 'calib_clip',[1 1], ...
    'noser_warmstart', struct('enable',false) );

%% ====== Lancement =======================================================
rec = simulate_inverse_slice(patient_this, slice_this, patient_prev, slice_prev, opts, iterDir);

%% ====== Sauvegarde finale ==============================================
save(recFile, '-struct', 'rec');

mis = NaN; if isfield(rec,'misfit') && ~isempty(rec.misfit), mis = rec.misfit; end
fprintf('\n[INVERSE-CONST] Reco terminée. Misfit relatif = %.6g\nFichier: %s\n', mis, recFile);
