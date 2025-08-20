% ============================
% File: Sensitivity_Sweep_Inverse.m
% ============================
% Balayage d'hyper-paramètres pour l'inverse EIT (OOEIT)
% Compare tv_alpha_scale, tv_ec, et (constNoise, relNoise)
% Sortie: tableau récap + 2 courbes (misfit vs alpha, misfit vs bruit)

clear; close all; clc;
addpath('src', genpath('src'));

% ---- cas utilisateur
patient_id = 's0011';
slice_this = 302;
slice_prev = 301;

% ---- grilles d'exploration (ajuste librement)
alphas = [1e-3, 3e-3, 1e-2, 3e-2, 1e-1];
ecs    = [0.10, 0.20, 0.30];
noises = [ % [constNoise, relNoise]
    5e-4  1e-2
    1e-3  3e-2
    2e-3  5e-2
];

% ---- autres options fixes
baseOpts = struct( ...
    'doWarp', true, ...
    'tv_beta', 1e-4, ...
    'clip', [0 inf], ...
    'savePlots', false, ...
    'dpi', 150 ...
);

% ---- chemins utiles pour évaluer les fonctionnels
rootThis = fullfile('Outputs', patient_id, sprintf('slice_%03d', slice_this));
packThis = fullfile(rootThis,'eit_pack.mat');
S = load(packThis);  % g,H,E,params,Imeas

% ---- collecte résultats
rows = {};
res  = [];
k = 0;

for ia = 1:numel(alphas)
    for ie = 1:numel(ecs)
        for in = 1:size(noises,1)

            k = k+1;
            opts = baseOpts;
            opts.tv_alpha_scale = alphas(ia);
            opts.tv_ec          = ecs(ie);
            opts.constNoise     = noises(in,1);
            opts.relNoise       = noises(in,2);

            % --- lance la reco (pas de PNG)
            rec = simulate_inverse_slice(patient_id, slice_this, slice_prev, opts);

            % --- évalue (en plus du misfit déjà présent) la valeur TV & data-term
            % reconstruire le modèle forward pour calculer data-term
            mesh = ForwardMesh1st(S.g*1e-3, S.H, S.E);
            fem  = EITFEM(mesh);
            fem.mode = 'current';
            fem.zeta = S.params.z_contact*ones(S.params.Ne,1);
            Iel = EITSim.buildTrigPattern(S.params.Ne, S.params.I_amp);
            fem.Iel = Iel(:);
            fem.Uel = S.Imeas(:);
            fem.SetInvGamma(opts.constNoise, opts.relNoise);

            Ufwd = fem.SolveForwardVec(rec.sigma_rec_nod);
            d    = Ufwd - S.Imeas(:);
            dataTerm = 0.5 * (d' * (fem.InvGamma * d));

            % valeur du prior TV avec mêmes hyper-paramètres
            tv = PriorTotalVariation(S.g, S.H, opts.tv_ec, 1, baseOpts.tv_beta);
            tv.alpha = tv.alpha * opts.tv_alpha_scale;     % IMPORTANT: même scaling
            tvVal   = tv.OptimizationFunction(rec.sigma_rec_nod);

            rows{k,1} = opts.tv_alpha_scale;
            rows{k,2} = opts.tv_ec;
            rows{k,3} = opts.constNoise;
            rows{k,4} = opts.relNoise;
            rows{k,5} = rec.data_misfit_rel;
            rows{k,6} = dataTerm;
            rows{k,7} = tvVal;
        end
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'tv_alpha_scale','tv_ec','constNoise','relNoise','misfit_rel','data_term','tv_value'});

% tri par misfit puis data_term
T = sortrows(T, {'misfit_rel','data_term'});

% ---- affichage + sauvegarde CSV
disp(T);
outCSV = fullfile(rootThis, 'sweep_inverse_results.csv');
writetable(T, outCSV);
fprintf('Résultats exportés: %s\n', outCSV);

% ---- quelques plots rapides
figure('Color','w','Name','Misfit vs alpha'); hold on;
for ie = 1:numel(ecs)
    idx = T.tv_ec == ecs(ie) & T.relNoise == noises(2,2) & T.constNoise == noises(2,1); % coupe à bruit "moyen"
    plot(T.tv_alpha_scale(idx), T.misfit_rel(idx), 'o-','DisplayName',sprintf('tv\\_ec=%.2f',ecs(ie)));
end
set(gca,'XScale','log'); grid on; xlabel('\alpha scale'); ylabel('misfit rel');
title('Misfit vs tv\_alpha\_scale (bruit moyen)'); legend('Location','best');

figure('Color','w','Name','Misfit vs bruit'); hold on;
labs = {};
for in = 1:size(noises,1)
    idx = T.constNoise==noises(in,1) & T.relNoise==noises(in,2) & T.tv_ec==ecs(2);
    plot(T.tv_alpha_scale(idx), T.misfit_rel(idx), 's-');
    labs{end+1} = sprintf('cN=%.0e, rN=%.0e',noises(in,1),noises(in,2)); %#ok<AGROW>
end
set(gca,'XScale','log'); grid on; xlabel('\alpha scale'); ylabel('misfit rel');
title('Misfit vs bruit (tv\_ec=0.20)'); legend(labs,'Location','best');
