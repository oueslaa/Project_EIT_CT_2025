% Driver_run_OOEIT_with_slices.m
% -------------------------------------------------------------------------
% Génère les slices OOEIT et affiche côte-à-côte avec les meshes DistMesh
% -------------------------------------------------------------------------

clear; close all; clc;

%% 1) SETTINGS
baseDir       = fileparts(mfilename('fullpath'));
addpath(fullfile('/Users/anis/Documents/StageInria/Code/','NIfTI_20140122'));

case_id       = 's0011';                              % <— Même case_id que test_tri2
Slices_folder = fullfile(baseDir,'Slices_folder');
if ~exist(Slices_folder,'dir'), mkdir(Slices_folder); end

% Répertoire des meshes triangulés
tri_dir = fullfile(baseDir,'triangulated_meshes',case_id);

% Charger le CT volume
ctData = load(fullfile(baseDir,'ct_mat.mat'));
ct_mat = ctData.ct_mat;

% Charger couleurs & noms complets
C           = load(fullfile(baseDir,'group_colors.mat'));
group_names = cellstr(C.group_names);  % ex. {'liver','colon',…,'soft'}
rgb_colors  = C.rgb_colors;
G           = numel(group_names);

% Paramètres de slice
z_min      = 201;
z_max      = 210;
slice_list = z_min:z_min;  % uniquement 201 ici

% Récupérer H×W depuis un masque fusionné (pour préallouer seg_img)
key0 = group_names{1};
if ~strcmpi(key0,'soft'), key0 = key0(1:4); end
nii0 = load_nii(fullfile(baseDir,'merged_masks',...
            sprintf('%s_mask_z%d_%d.nii.gz', key0, z_min, z_max)));
[h, w, ~] = size(nii0.img);


%% 2) BOUCLE SUR CHAQUE SLICE
for iz = slice_list
    fprintf('--- Slice %d ---\n', iz);

    %% A) CALCUL DU MAILLAGE OOEIT
    slice_raw  = double(ct_mat.img(:,:,iz));
    slice_norm = (slice_raw - min(slice_raw(:))) / range(slice_raw(:));

    BW   = slice_raw > -500;
    L    = bwlabel(BW,8);
    cnts = histcounts(L(L>0),1:(max(L(:))+1));
    [~, idx] = max(cnts);
    mask = imfill(L==idx,'holes');

    down = 0.25;
    mask_q  = imresize(mask,down,'nearest');
    slice_q = imresize(slice_norm,down);
    [nx, ny] = size(mask_q);

    [jj, kk] = find(mask_q);
    shp       = alphaShape(kk,jj);
    facets0   = shp.boundaryFacets;
    pts0      = shp.Points;

    ang = atan2(pts0(facets0(:,1),2)-nx/2, pts0(facets0(:,1),1)-ny/2);
    fitx = fit(ang, pts0(facets0(:,1),1)-ny/2, "fourier8");
    fity = fit(ang, pts0(facets0(:,1),2)-nx/2, "fourier8");
    pv   = [1.1*fitx(ang)+ny/2, 1.1*fity(ang)+nx/2];

    fd = {'l_dpolygon',[],pv};
    fh = @(p) ones(size(p,1),1);

    [p, t] = distmesh(fd, fh, 2, [0,0;nx,ny], pv);

    TR      = triangulation(t,p);
    edges   = TR.edges;
    attach  = TR.edgeAttachments(edges);
    bd_e    = edges(cellfun(@(c)numel(c)==1,attach),:);

    Vq  = interp2(1:ny,1:nx,slice_q,p(:,1),p(:,2),'nearest');
    sig = Vq;
    for ia = 1:size(bd_e,1)
        i0 = bd_e(ia,1);
        d  = hypot(p(:,1)-p(i0,1), p(:,2)-p(i0,2));
        sig(d<=3) = 0.5;
    end

    %% B) RECONSTRUCTION SEGMENTATION RGB
    seg_img = 255 * ones(h, w, 3, 'uint8');
    for k = 1:G
        key = group_names{k};
        if strcmpi(key,'soft')
            shortkey = 'soft';
        else
            shortkey = key(1:4);
        end

        nii = load_nii(fullfile(baseDir,'merged_masks',...
               sprintf('%s_mask_z%d_%d.nii.gz', shortkey, z_min, z_max)));
        sl     = iz - z_min + 1;
        mask2d = nii.img(:,:,sl)>0;

        for c = 1:3
            layer = seg_img(:,:,c);
            layer(mask2d) = uint8(255 * rgb_colors(k,c));
            seg_img(:,:,c) = layer;
        end
    end

    %% C) AFFICHAGE CÔTE-À-CÔTE
    fig = figure('Units','normalized','Position',[0.1 0.1 0.8 0.6]);
    tl  = tiledlayout(1,2, ...
         'TileSpacing','compact', ...
         'Padding',    'compact');

    % 1) OOEIT mesh
    ax1 = nexttile(tl,1);
    trimesh(TR.ConnectivityList, p(:,1), p(:,2), sig, ...
            'EdgeColor','none', 'FaceColor','interp', 'FaceAlpha',0.6);
    axis(ax1,'equal','off','ij');
    title(ax1, sprintf('OOEIT Slice %d', iz));

    % 2) Segmentation + DistMesh
    ax2 = nexttile(tl,2);
    imshow(seg_img,'Parent',ax2);
    axis(ax2,'image','off'); hold(ax2,'on');
    title(ax2, sprintf('Segmentation + Meshs – Slice %d', iz));

    % Superpose les maillages
    axes(ax2);  % on se place sur ax2
    for k = 1:G
        key = group_names{k};
        if strcmpi(key,'soft'), shortkey='soft';
        else shortkey = key(1:4);
        end
        mesh_mat = fullfile(tri_dir, sprintf('slice_%d',iz), shortkey,'mesh.mat');
        if exist(mesh_mat,'file')
            S = load(mesh_mat);
            triplot(S.t, S.p(:,1), S.p(:,2), 'k', 'LineWidth',0.5);
        end
    end
    hold(ax2,'off');

    % Légende
    ph = gobjects(G,1);
    for k = 1:G
        ph(k) = patch(ax2, nan, nan, rgb_colors(k,:), ...
                       'EdgeColor','none', 'DisplayName',group_names{k});
    end
    legend(ax2, ph, group_names, ...
           'Location','southoutside', ...
           'Orientation','horizontal', ...
           'NumColumns',ceil(G/5), ...
           'Box','off');

    drawnow; pause(0.1);

    %% D) SAUVEGARDE OOEIT
    save(fullfile(Slices_folder, sprintf('slice_OOEIT_%d.mat', iz)), ...
         'p','t','sig');
    fprintf('→ Saved slice_OOEIT_%d.mat\n', iz);
end

fprintf('✅ Toutes les slices sont dans %s\n', Slices_folder);
