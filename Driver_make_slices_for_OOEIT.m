% Driver_make_slices_for_OOEIT.m
% -------------------------------------------------------------------------
% Génère des tranches 2D à partir de ct_mat.img, construit g, H, elfaces, sig,
% puis les enregistre sous forme de slice_OOEIT_<z>.mat dans Slices_folder.
% Affiche en plus, pour chaque slice, la carte OOEIT à gauche et l'image
% segmentée colorée à droite, avec légende des groupes.
% -------------------------------------------------------------------------

%% 1) Réglages de base
clear; close all; clc;

% Localisation du script et ajout de DistMesh
baseDir     = fileparts(mfilename('fullpath'));
distmeshDir = fullfile(baseDir, '..', 'distmesh-master');
if ~exist(distmeshDir, 'dir')
    error('Le dossier DistMesh introuvable : %s', distmeshDir);
end
addpath(genpath(distmeshDir));

% Chemin vers le dossier des slices
Slices_folder = fullfile(baseDir, 'Slices_folder');
if ~exist(Slices_folder, 'dir')
    mkdir(Slices_folder);
    fprintf('Création du dossier : %s\n', Slices_folder);
end

% Chargement du volume CT
ctData = load(fullfile(baseDir, 'ct_mat.mat'));
ct_mat = ctData.ct_mat;

%% 2) Paramètres de découpe
ns = 201;   % slice de début
nf = 210;   % slice de fin

% — charger les légendes de segmentation (générées depuis Python)
c = load(fullfile(baseDir, 'group_colors.mat'));
group_names = c.group_names;     % peut être char array ou cell array
rgb_colors  = c.rgb_colors;      % matrice N×3

% Convertir group_names en cell array si besoin
if ischar(group_names) || isstring(group_names)
    group_names = cellstr(group_names);
end

% Vérification / ajustement de la taille
G1 = numel(group_names);
G2 = size(rgb_colors,1);
if G1 ~= G2
    warning('Group names (%d) et RGB rows (%d) diffèrent : on tronque.', G1, G2);
    G = min(G1, G2);
    group_names = group_names(1:G);
    rgb_colors  = rgb_colors(1:G, :);
else
    G = G1;
end

%% 3) Boucle de traitement
for iz = ns:nf
    fprintf('--- Traitement slice %d ---\n', iz);

    % Extraction et normalisation de l'image
    slice_raw  = double(ct_mat.img(:,:,iz));
    minI = min(slice_raw(:));
    slice_pos  = slice_raw - minI;
    slice_norm = slice_pos / max(slice_pos(:));

    % Segmentation du blob principal
    BW     = slice_raw > -500;
    L      = bwlabel(BW, 8);
    counts = histcounts(L(L>0), 1:(max(L(:))+1));
    [~,idx]= max(counts);
    mask   = imfill(L==idx,'holes');
    mask_q = imresize(mask,0.25,'nearest');
    slice_q= imresize(slice_norm,0.25);

    % Construction de l'alphaShape pour la frontière
    [jj,kk] = find(mask_q);
    shp     = alphaShape(kk,jj);
    facets0 = shp.boundaryFacets;
    pts0    = shp.Points;

    % Lissage Fourier
    nx = size(mask_q,1); ny = size(mask_q,2);
    ii = facets0(:,1);
    ang = atan2(pts0(ii,2)-nx/2, pts0(ii,1)-ny/2);
    fitx = fit(ang, pts0(ii,1)-nx/2, "fourier8");
    fity = fit(ang, pts0(ii,2)-ny/2, "fourier8");
    fx_smooth = 1.1*fitx(ang) + nx/2;
    fy_smooth = 1.1*fity(ang) + ny/2;
    pv = [fx_smooth, fy_smooth];

    % DistMesh
    fd = { 'l_dpolygon', [], pv };
    fh = @(p) ones(size(p,1),1);
    [p,t] = distmesh(fd, fh, 2, [0,0; nx,ny], pv);

    % Triangulation
    TR = triangulation(t,p);
    edges = TR.edges;
    attach = TR.edgeAttachments(edges);
    bd_edges = edges(cellfun(@(c)numel(c)==1,attach),:);

    % Ordonnancement frontière
    facets_ord = bd_edges(1,:);
    bd = bd_edges;
    for k=2:size(bd,1)
        prev = facets_ord(end,2);
        idx1 = find(bd(:,1)==prev,1);
        idx2 = find(bd(:,2)==prev,1);
        if ~isempty(idx1)
            next = bd(idx1,:);
            bd(idx1,:) = nan;
        else
            next = bd(idx2,[2 1]);
            bd(idx2,:) = nan;
        end
        facets_ord(end+1,:) = next;
    end

    % Découpage électrodes
    nel = 16;
    nfa = size(facets_ord,1);
    vec = floor(linspace(1,nfa,nel+1));
    elfaces = cell(nel,1);
    for ie=1:nel
        startIdx = vec(ie);
        endIdx   = min(vec(ie)+3, nfa);
        elfaces{ie} = facets_ord(startIdx:endIdx,:);
    end

    % Interpolation sigma
    xvec = linspace(1,78,size(slice_q,1));
    yvec = linspace(1,78,size(slice_q,2));
    Vq = interp2(xvec,yvec,slice_q,p(:,1),p(:,2),'nearest');

    % Masquage bords
    sig = Vq;
    bdy_gap = 3;
    for ia=1:size(facets_ord,1)
        d = hypot(p(:,1)-p(facets_ord(ia,1),1), p(:,2)-p(facets_ord(ia,1),2));
        sig(d<=bdy_gap) = 0.5;
    end

    g = p;
    H = TR.ConnectivityList;

    %% Affichage côte-à-côte
    fig = figure('Units','normalized','Position',[0.1 0.1 0.8 0.6]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % 1) Carte OOEIT
    ax1 = nexttile;
    h = trimesh(H(:,1:3),g(:,1),g(:,2),sig);
    set(h,'FaceColor','interp','FaceAlpha',0.6,'EdgeColor','none');
    axis(ax1,'equal','off');
    view(ax1,2);
    title(ax1,sprintf('OOEIT Slice %d',iz));
    colorbar(ax1);

    % 2) Segmentation colorée (miroir horizontal)
    ax2 = nexttile;
    img_seg_orig = imread(fullfile(baseDir, 'segmentation_slices', sprintf('seg_%d.png', iz)));
    img_seg      = flipud(img_seg_orig);        % flip haut↔bas (miroir horizontal)
    imshow(img_seg, 'Parent', ax2);
    axis(ax2, 'off');
    title(ax2, sprintf('Segmentation colorée – Slice %d', iz));


    

    % 3) Légende sous la segmentation
    hold(ax2,'on');
    pHandle = gobjects(G,1);
    for k=1:G
        pHandle(k) = patch(ax2,nan,nan,rgb_colors(k,:),'EdgeColor','none','DisplayName',group_names{k});
    end
    hold(ax2,'off');
    leg = legend(ax2,pHandle,group_names,...
                 'Location','southoutside',...
                 'Orientation','horizontal',...
                 'NumColumns',ceil(G/5),...
                 'Box','off');
    % Forcer la légende sous les deux tuiles (si supported)
    try
        leg.Layout.Tile = 'south';
    end

    drawnow;
    pause(0.2);

    %% Sauvegarde des données OOEIT
    sliceName = sprintf('slice_OOEIT_%d.mat', iz);
    save(fullfile(Slices_folder,sliceName),'g','H','elfaces','sig');
    fprintf('→ Enregistré : %s\n', sliceName);
end

fprintf('Terminé : toutes les slices sont dans %s\n',Slices_folder);
