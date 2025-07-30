%% compare_slices_eit.m
clear; close all; clc; format compact;
tstart = tic;

MeshDataDir = '/Users/anis/Documents/StageInria/Code/Project_EIT_CT_2025-main/MeshData';
zvec = 201:205;              % Slices à comparer
iz_sol = 202;                % Référence
Ne = 16;                     % Nombre d'électrodes
nbest = 2;                   % Nombre de meilleurs slices à afficher

% ----------- Chargement du mesh de référence et calcul du signal ----------
S = load(fullfile(MeshDataDir, sprintf('mesh_slice%03d.mat', iz_sol)), ...
          'g','H','elfaces','sig');
g = S.g; H = S.H; elfaces = S.elfaces; sig = S.sig;

g = g / 1000;   % Conversion mm -> m pour le solveur EIT

if isempty(sig) || size(sig,1) ~= size(g,1)
    error('Conductivité absente ou non interpolée pour la tranche de référence !');
end

mesh = ForwardMesh1st(g, H, elfaces);
solver = EITFEM(mesh);
solver.mode = 'potential';
Imeas = solver.SolveForwardVec(sig);

signal_sol = Imeas;

% Suppression des mesures négatives et de leurs voisins
for n = 0:Ne-1
    idx = (1:Ne) + n*Ne;
    ii = find(signal_sol(idx) < 0);
    iip = mod(ii, Ne) + 1;
    iim = mod(ii-2, Ne) + 1;
    signal_sol(idx(ii)) = NaN;
    signal_sol(idx(iip)) = NaN;
    signal_sol(idx(iim)) = NaN;
end

signal = cell(1, max(zvec));
nodes  = cell(1, max(zvec));
elem   = cell(1, max(zvec));
cond   = cell(1, max(zvec));

% ----------- Boucle sur toutes les slices à comparer ----------
for iz = zvec
    S = load(fullfile(MeshDataDir, sprintf('mesh_slice%03d.mat', iz)), ...
              'g','H','elfaces','sig');
    g = S.g; H = S.H; elfaces = S.elfaces; sig = S.sig;

    g = g / 1000;   % Conversion mm -> m pour le solveur EIT

    if isempty(sig) || size(sig,1) ~= size(g,1)
        warning('Slice %d : sig absent ou incompatible, ignorée', iz);
        continue;
    end
    mesh = ForwardMesh1st(g, H, elfaces);
    solver = EITFEM(mesh);
    solver.mode = 'potential';
    Imeas = solver.SolveForwardVec(sig);
    signal{iz} = Imeas;
    nodes{iz}  = g;
    elem{iz}   = H;
    cond{iz}   = sig;
end

% ----------- Calcul de l'erreur L2 ----------
err_vec = zeros(numel(zvec),1);
innan = ~isnan(signal_sol);

for i = 1:numel(zvec)
    iz = zvec(i);
    Imeas = signal{iz};
    if isempty(Imeas)
        err_vec(i) = NaN; continue;
    end
    for n = 0:Ne-1
        idx = (1:Ne) + n*Ne;
        ii = find(Imeas(idx) < 0);
        iip = mod(ii, Ne) + 1;
        iim = mod(ii-2, Ne) + 1;
        Imeas(idx(ii)) = NaN;
        Imeas(idx(iip)) = NaN;
        Imeas(idx(iim)) = NaN;
    end
    err_vec(i) = norm(Imeas(innan) - signal_sol(innan));
end

% ----------- Affichage du graphe de l'erreur ----------
figure;
plot(zvec, err_vec, 'x-', 'LineWidth',2);
xlabel('Indice de slice');
ylabel('Erreur L2');
title('Erreur L2 des tranches comparées');
grid on;

% ----------- Affichage des 2 meilleures slices ----------
[~, best_idx] = mink(err_vec, nbest);
for k = 1:nbest
    iz_best = zvec(best_idx(k));
    g = nodes{iz_best};
    H = elem{iz_best};
    Imeas = signal{iz_best};
    sig = cond{iz_best};

    figure('Position',[100 100 1800 600]);
    subplot(1,3,1);
    plot(Imeas,'-o');
    xlabel('Mesure #');
    ylabel('Signal');
    title(sprintf('Tranche %d — signal simulé', iz_best));
    grid on;

    subplot(1,3,2);
    patch('Faces',H,'Vertices',g, ...
        'FaceVertexCData',sig, ...
        'FaceColor','interp','EdgeAlpha',0.2,'FaceAlpha',0.95);
    axis equal tight;
    xlabel('X (m)'); ylabel('Y (m)');
    colorbar; colormap parula;
    title(sprintf('Mesh & CT interpolée — slice %d', iz_best));

    % Affichage PNG preview à droite (si dispo)
    subplot(1,3,3);
    try
        img_png = imread(fullfile(MeshDataDir, sprintf('preview_slice%03d.png', iz_best)));
        imshow(img_png); axis off;
        title('CT + Contours + Mesh');
    catch
        text(0.2,0.5,'Preview PNG not found');
        axis off;
    end
end

fprintf('Comparaison terminée en %.2f secondes.\n', toc(tstart));
