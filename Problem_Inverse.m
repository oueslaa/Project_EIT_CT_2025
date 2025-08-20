% ============================
% File: Problem_Inverse.m
% ============================
% Reconstruction EIT de slice_this en partant de slice_prev.
% Pipeline robuste : multi-stages GN (WLS + TV + L2-to-init) + soft-snap.
% Rien à modifier dans OOEIT, le forward, ni init_from_prev_slice_discrete.

clear; close all; clc;
addpath('src', genpath('src'));
addpath(genpath('/Users/anis/Documents/StageInria/Code/OOEIT-main'));  % adapte si besoin

% Afficher les figures (plots finaux)
set(groot,'defaultFigureVisible','on');

% -------- Choix patient / slices --------
patient_id  = 's0011';
slice_this  = 301;   % cible à reconstruire
slice_prev  = 302;   % source pour l'initialisation

% --------- Paramètres simples (reco robuste une passe) ----------
opts = struct();

% Bornes physiques
opts.clip = [0.05, 0.50];

% Init (discrète + léger nettoyage) + recalib d'amplitude
opts.init = struct();
opts.init.discrete       = true;
opts.init.calib_enable   = true;
opts.init.calib_clip     = [0.95, 1.05];           % un peu plus serré
opts.init.quant          = struct('centers',[0.05 0.15 0.25 0.45], 'clean_iters',2);
opts.doWarp              = true;

% Pondération bruit (WLS) — remets les défauts "peu de bruit"
opts.noise = struct('const', 2e-5, 'rel', 2e-2);

% Régularisation — celles qui marchaient bien
opts.tv_ec    = 0.20;
opts.tv_beta  = 7e-5;
opts.tv_alpha = 1.8e-4;      % TV un poil plus forte que 2e-4, moins que 4e-4
opts.l2_alpha = 1e-2;      % L2-to-init assez ferme (aide à garder l'os à 0.05)

% Gauss-Newton (une passe)
opts.maxIter        = 8;   % 5–8 OK
opts.maxIterInLine  = 25;

% Pas de snap final
opts.snap_final_to_levels = false;
opts.min_patch_tri = 0;

% Sorties
opts.savePlots = true;
opts.dpi       = 300;

% -------- Run --------
rec = simulate_inverse_slice(patient_id, slice_this, slice_prev, opts);

fprintf('\nReco slice %d terminée.\n  Misfit relatif = %.4g\n  Résultats : %s\n', ...
        slice_this, rec.data_misfit_rel, rec.paths.outMat);
