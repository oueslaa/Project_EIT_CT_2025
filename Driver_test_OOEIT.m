%% 1) Réglages de base
clear; close all; clc;

baseDir = fileparts(mfilename('fullpath'));

ooeitDir = fullfile(baseDir, '..', 'OOEIT-main');
if ~exist(ooeitDir, 'dir')
    error('Le dossier OOEIT introuvable : %s', ooeitDir);
end
addpath(genpath(ooeitDir));

libDir  = fullfile(baseDir, 'meshes_mat'); % dossier où tu as sauvegardé les mesh_*.mat
refFile = 's0011_mesh.mat';                % fichier de référence (g,H,sigma)

b = [0.3, 0.15, 0.01, 0.2, 0.6, 0.8, 0.1, 0.25, 0.05];
n_el = 16;

%% 2) Charger la référence
Sref = load(fullfile(libDir,refFile));
g_ref = Sref.g; H_ref = Sref.H + 1; sig_ref = Sref.sigma;

c = zeros(size(sig_ref));
for i=1:numel(b), c(sig_ref==i) = b(i); end
sig_ref = c;

E_ref = cell(n_el,1);
for L=1:n_el
    E_ref{L} = Sref.(sprintf('E%d',L)) + 1;
end

fmesh_ref  = ForwardMesh1st(g_ref,H_ref,E_ref);
solver_ref = EITFEM(fmesh_ref); solver_ref.mode='current';
Imeas_ref  = solver_ref.SolveForwardVec(sig_ref);

%% 3) Boucle sur les fichiers de maillage
files = dir(fullfile(libDir,'mesh_*.mat'));
nF    = numel(files);

results = table('Size',[nF 4], ...
                'VariableTypes',{'string','double','double','double'}, ...
                'VariableNames',{'File','NbTriangles','ForwardTime','ErrRef'});

for k=1:nF
    fname = files(k).name;
    Smesh = load(fullfile(libDir,fname));

    % ⚠️ Charger le maillage PDE Toolbox
    g = Smesh.mesh.Nodes';       % Transposé
    H = Smesh.mesh.Elements';    % Transposé
    sig = ones(size(H,1),1);     % conductivité uniforme (par défaut)

    % Placer 16 électrodes fictives
    theta = linspace(0,2*pi,n_el+1);
    x_el  = mean(g(:,1)) + 150*cos(theta(1:end-1))';
    y_el  = mean(g(:,2)) + 150*sin(theta(1:end-1))';
    E = cell(n_el,1);
    for L=1:n_el
        [~,idx] = min(vecnorm(g-[x_el(L),y_el(L)],2,2));
        E{L} = idx;
    end

    % Forward solve avec mesure du temps
    fmesh  = ForwardMesh1st(g,H,E);
    solver = EITFEM(fmesh); solver.mode='current';

    tic;
    Imeas = solver.SolveForwardVec(sig);
    tSolve = toc;

    errRef = norm(Imeas - Imeas_ref);

    results(k,:) = {fname,size(H,1),tSolve,errRef};

    fprintf('%2d/%2d %-25s | Tri=%5d | Time=%.3fs | ErrRef=%.3e\n', ...
        k,nF,fname,size(H,1),tSolve,errRef);
end

%% 4) Affichage des résultats triés par temps
results = sortrows(results,'ForwardTime');
disp(results);

figure;
bar(results.ForwardTime);
xticklabels(results.File);
xtickangle(45);
ylabel('Temps forward solve (s)');
title('Comparaison des maillages PDE Toolbox');
grid on;
