%% mesh_extract_slice_colored_triangles.m
clear; close all; clc;

patient_id = 's0011';
z_slice    = 301;
meshFile   = fullfile('MeshData', sprintf('mesh_slice%03d.mat', z_slice));
if ~exist(meshFile, 'file')
    error('Fichier de maillage introuvable : %s', meshFile);
end

% Charger le mesh brut
load(meshFile, 'g', 'H');

% Récupérer polygones
warnState = warning('off','MATLAB:polyshape:checkAndSimplify');
[domain, shapes, ext_smooth, info] = extract_slice_polys(patient_id, z_slice);
warning(warnState);

types = {'SoftTissue','Heart','Lung','Trachea','Bone','Other'};
nG    = numel(types);

% Centroides des triangles
centroids = ( g(H(:,1),:)+g(H(:,2),:)+g(H(:,3),:) )/3;

% Attribution groupe
M = size(H,1);
triGroup = zeros(M,1);
for t = 1:M
    pt = centroids(t,:);
    assigned = false;
    for k = 2:nG
        if isinterior(shapes.(types{k}), pt(1), pt(2))
            triGroup(t) = k;
            assigned = true;
            break;
        end
    end
    if ~assigned && isinterior(domain, pt(1), pt(2))
        triGroup(t) = 1; % soft tissue
    end
end

% Sauvegarder complet
save(meshFile, 'g', 'H', 'triGroup', 'shapes', 'domain', '-append');
fprintf('Mesh %s mis à jour avec triGroup, shapes, domain.\n', meshFile);
