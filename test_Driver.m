% Driver_make_slices_for_OOEIT.m
% -------------------------------------------------------------------------
% Génère des tranches OOEIT à partir de ct_mat.mat,
% reconstruit l'image segmentée slice par slice depuis les masques NIfTI fusionnés,
% et affiche côte à côte :
%   • à gauche : votre carte OOEIT
%   • à droite : segmentation colorée + meshes DistMesh
% -------------------------------------------------------------------------

clear; close all; clc;

%% 1) Réglages de base
baseDir       = fileparts(mfilename('fullpath'));
addpath(fullfile(baseDir, 'NIfTI_20140122'));  % pour load_nii
Slices_folder = fullfile(baseDir, 'Slices_folder');
if ~exist(Slices_folder,'dir'), mkdir(Slices_folder); end

% Charger le volume CT
ctData = load(fullfile(baseDir, 'ct_mat.mat'));
ct_mat = ctData.ct_mat;

% Charger couleurs & noms de groupes
C           = load(fullfile(baseDir, 'group_colors.mat'));
group_names = cellstr(C.group_names);
rgb_colors  = C.rgb_colors;
G           = numel(group_names);

% Paramètres de slice
z_min      = 201;
z_max      = 210;
slice_list = z_min:z_min;

% Répertoire des meshes produits par test_tri2.m
tri_dir = fullfile(baseDir, 'triangulated_meshes', 'Project_EIT_CT_2025-main');

% Préparer h,w depuis un masque NIfTI
nii0 = load_nii(fullfile(baseDir,'merged_masks', ...
         sprintf('%s_mask_z%d_%d.nii.gz', group_names{1}, z_min, z_max)));
[h, w, ~] = size(nii0.img);

%% 2) Boucle sur chaque slice
for iz = slice_list
    fprintf('--- Traitement slice %d ---\n', iz);

    % A) Calcul OOEIT
    slice_raw  = double(ct_mat.img(:,:,iz));
    slice_norm = (slice_raw - min(slice_raw(:))) / range(slice_raw(:));
    BW     = slice_raw > -500;
    L      = bwlabel(BW,8);
    counts = histcounts(L(L>0), 1:(max(L(:))+1));
    [~,idx]= max(counts);
    mask   = imfill(L==idx,'holes');

    downscale = 0.25;
    mask_q    = imresize(mask,   downscale,'nearest');
    slice_q   = imresize(slice_norm, downscale);
    nx = size(mask_q,1);
    ny = size(mask_q,2);

    [jj,kk] = find(mask_q);
    shp     = alphaShape(kk,jj);
    facets0 = shp.boundaryFacets;
    pts0    = shp.Points;

    ang  = atan2( pts0(facets0(:,1),2)-nx/2, pts0(facets0(:,1),1)-ny/2 );
    fitx = fit(ang, pts0(facets0(:,1),1)-ny/2, "fourier8");
    fity = fit(ang, pts0(facets0(:,1),2)-nx/2, "fourier8");
    pv   = [1.1*fitx(ang)+ny/2, 1.1*fity(ang)+nx/2];

    fd = {'l_dpolygon',[],pv};
    fh = @(p) ones(size(p,1),1);
    [p,t] = distmesh(fd, fh, 2, [0,0; nx,ny], pv);

    TR       = triangulation(t,p);
    edges    = TR.edges;
    attach   = TR.edgeAttachments(edges);
    bd_edges = edges(cellfun(@(c)numel(c)==1,attach),:);

    Vq  = interp2(1:ny, 1:nx, slice_q, p(:,1), p(:,2),'nearest');
    sig = Vq;
    for ia = 1:size(bd_edges,1)
        idxb = bd_edges(ia,1);
        d    = hypot(p(:,1)-p(idxb,1), p(:,2)-p(idxb,2));
        sig(d<=3) = 0.5;
    end

    %% Reconstruction de l'image segmentée RGB
    seg_img = zeros(h, w, 3, 'uint8');
    for k = 1:G
        key = group_names{k};
        nii = load_nii(fullfile(baseDir,'merged_masks', ...
               sprintf('%s_mask_z%d_%d.nii.gz', key, z_min, z_max)));
        sl  = iz - z_min + 1;
        mask2d = nii.img(:,:,sl) > 0;
        for c = 1:3
            layer = seg_img(:,:,c);
            layer(mask2d) = uint8(255 * rgb_colors(k,c));
            seg_img(:,:,c) = layer;
        end
    end
    mask_any = any(seg_img,3);
    seg_img(repmat(~mask_any,[1 1 3])) = 255;  % fond blanc

    %% Affichage côte-à-côte
    fig = figure('Units','normalized','Position',[0.1 0.1 0.8 0.6]);
    tl  = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % 1) OOEIT
    ax1 = nexttile(tl,1);
    hMesh = trimesh(TR.ConnectivityList, p(:,1), p(:,2), sig, 'Parent', ax1);
    set(hMesh,'FaceColor','interp','FaceAlpha',0.6,'EdgeColor','none');
    axis(ax1,'equal','off','ij');
    view(ax1,2);
    title(ax1, sprintf('OOEIT Slice %d', iz));

    % 2) Segmentation + Meshs
    ax2 = nexttile(tl,2);
    imshow(seg_img,'Parent',ax2);
    axis(ax2,'image','off'); hold(ax2,'on');
    title(ax2, sprintf('Segmentation + Meshs – Slice %d', iz));

    % Superposition des maillages DistMesh en noir
    for k = 1:G
        mesh_mat = fullfile(tri_dir, sprintf('slice_%d',iz), group_names{k}, 'mesh.mat');
        if exist(mesh_mat,'file')
            S = load(mesh_mat);
            triplot(S.t, S.p(:,1), S.p(:,2), 'Parent',ax2, 'Color','k','LineWidth',0.5);
        end
    end

    % Légende avec les vraies couleurs
    pHandle = gobjects(G,1);
    for k=1:G
        pHandle(k) = patch(ax2, nan, nan, rgb_colors(k,:), ...
                           'EdgeColor','none', 'DisplayName', group_names{k});
    end
    legend(ax2, pHandle, group_names, ...
           'Location','southoutside','Orientation','horizontal', ...
           'NumColumns',ceil(G/5),'Box','off');
    hold(ax2,'off');

    drawnow; pause(0.1);

    %% Sauvegarde OOEIT
    sliceName = sprintf('slice_OOEIT_%d.mat', iz);
    save(fullfile(Slices_folder, sliceName), 'p','t','sig');
    fprintf('→ Enregistré : %s\n', sliceName);
end

fprintf('✅ Toutes les slices sont dans %s\n', Slices_folder);
