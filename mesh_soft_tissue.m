% Script: mesh_extract_slice_polys_mask2mesh.m
% Génère un maillage 2D à partir des contours extraits par extract_slice_polys
% et colore chaque groupe d'organes par la même couleur

clear; close all; clc;

% Paramètres utilisateur
targetSize  = 5;    % Taille des éléments de maillage en mm
minArea_mm2 = 300;  % Seuil minimal d'aire en mm²
patient_id  = 's0011';
z_slice     = 301;
outDir      = 'MeshData';
if ~exist(outDir, 'dir'), mkdir(outDir); end

tol = 1e-6; % tolérance points trop proches

% 1) Extraction des contours via extract_slice_polys
warnState = warning('off','MATLAB:polyshape:checkAndSimplify');
[~, shapes, ext_smooth, info] = extract_slice_polys(patient_id, z_slice);
warning(warnState);

% 2) Collecte et nettoyage des contours
cleanContours = {};
names        = {};

types = {'Ext','Heart','Lung','Trachea','Bone','Other'};  % ordre des groupes

% 2.1) Contour extérieur du soft tissue
P_ext = ext_smooth;  % [Nx2]
dP    = [inf; hypot(diff(P_ext(:,1)), diff(P_ext(:,2)))];
P_ext(dP < tol, :) = [];
P_ext = unique(P_ext, 'rows', 'stable');
if ~isequal(P_ext(1,:), P_ext(end,:)), P_ext(end+1,:) = P_ext(1,:); end
P_ext(end,:) = [];
if size(P_ext,1) < 3, error('Contour extérieur insuffisant.'); end
cleanContours{end+1} = P_ext;
names{end+1} = 'Ext';

% 2.2) Contours des organes
orgFields = fieldnames(shapes);
for k = 1:numel(orgFields)
    shp = shapes.(orgFields{k});
    if shp.NumRegions == 0, continue; end
    regs = regions(shp);
    grp = orgFields{k};
    for r = 1:numel(regs)
        P = regs(r).Vertices;
        dP = [inf; hypot(diff(P(:,1)), diff(P(:,2)))];
        P(dP < tol, :) = [];
        P = unique(P, 'rows', 'stable');
        if ~isequal(P(1,:), P(end,:)), P(end+1,:) = P(1,:); end
        P(end,:) = [];
        if size(P,1) < 3 || polyarea(P(:,1),P(:,2)) < minArea_mm2, continue; end
        cleanContours{end+1} = P;
        names{end+1}        = grp;
    end
end
if isempty(cleanContours), error('Aucun contour valide.'); end

% 3) Construction des matrices gd, sf, ns
nS = numel(cleanContours);
vk = cellfun(@(C) size(C,1), cleanContours);
maxk = max(vk);
gd = zeros(2 + 2*maxk, nS);
for i = 1:nS
    P = cleanContours{i}; k = size(P,1);
    gd(1,i)          = 2;
    gd(2,i)          = k;
    gd(3:2+k, i)     = P(:,1);
    gd(3+k:2+2*k, i) = P(:,2);
end
sf = strjoin(arrayfun(@(i) sprintf('P%d', i), 1:nS, 'UniformOutput', false), '+');
ns = char(arrayfun(@(i) sprintf('P%d', i), 1:nS, 'UniformOutput', false))';

% 4) Création du modèle PDE et géométrie
disp('Création géométrie...');
model = createpde();
try
    [dl, ~] = decsg(gd, sf, ns);
    geometryFromEdges(model, dl);
catch ME
    error('decsg échoué : %s', ME.message);
end

% 5) Génération du maillage
disp('Génération du maillage...');
try
    msh = generateMesh(model, ...
        'Hmax',           targetSize, ...
        'Hmin',           targetSize, ...
        'Hgrad',          1,           ...
        'GeometricOrder', 'linear',    ...
        'MesherVersion',  'preR2013a', ...
        'Jiggle',         {'on','minimum'});
catch ME
    error('generateMesh échoué : %s', ME.message);
end

g = msh.Nodes';    % [N x 2]
H = msh.Elements'; % [M x 3]

% 6) Affichage du résultat par groupe
figure; hold on; axis equal;
% tracer maillage
triplot(H, g(:,1), g(:,2), 'Color', [0 0 0 0.4]);
% préparer couleurs pour chaque type
groups = types;
cols   = lines(numel(groups));
% tracer contours par groupe
for i = 1:nS
    P = cleanContours{i};
    grp = names{i};
    idx = find(strcmp(groups, grp));
    plot(P(:,1), P(:,2), 'LineWidth', 1.5, 'Color', cols(idx, :));
end
legendEntries = [{'Mesh'}, groups];
legend(legendEntries, 'Location', 'BestOutside');
title(sprintf('Maillage slice %d avec contours groupés', z_slice));

% 7) Sauvegarde des données
outFile = fullfile(outDir, sprintf('mesh_slice%03d.mat', z_slice));
save(outFile, 'g', 'H');
fprintf('Maillage enregistré dans %s\n', outFile);
