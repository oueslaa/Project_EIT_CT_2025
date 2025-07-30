%% mesher.m : Génère MeshData/mesh_sliceXXX.mat à partir du masque binaire NIfTI + CT
clear; close all; clc;
tic_all = tic;

% ===== PARAMÈTRES UTILISATEUR =====
nii_path = '/Users/anis/Documents/StageInria/Code/Project_EIT_CT_2025-main/merged_masks/soft_mask_z201_210.nii.gz';
ct_path  = '/Users/anis/Documents/StageInria/Code/Project_EIT_CT_2025-main/Data_set/s0011/ct.nii.gz';
outDir   = 'MeshData';
if ~exist(outDir, 'dir'), mkdir(outDir); end
K_keep      = 12;        % lissage FFT contours
downsample  = 2;
minArea_mm2 = 300;       % aire min contour (mm²)
baseSlice   = 201;       % index du premier slice
targetSize  = 5;         % taille cible triangle (mm)
Ne          = 16;        % nombre d'électrodes

% ====== Lecture volumes ======
info  = niftiinfo(nii_path);
vol3d = niftiread(info);    % [ny,nx,nz]
dx = info.PixelDimensions(1);
dy = info.PixelDimensions(2);
[ny,nx,nz] = size(vol3d);

info_ct  = niftiinfo(ct_path);
vol_ct   = niftiread(info_ct);

fprintf('Segmentation: %dx%dx%d (%s)\n', ny,nx,nz, nii_path);
fprintf('CT: %dx%dx%d (%s)\n', size(vol_ct), ct_path);

for k = 1:nz
    t0 = tic;
    sliceNum = baseSlice + k - 1;
    mask2D   = squeeze(vol3d(:,:,k)) > 0;

    % === Extraction et nettoyage des contours ===
    [B,~] = bwboundaries(mask2D, 'holes');
    cleanContours = {};
    for i = 1:numel(B)
        pts = B{i}(:,[2,1]).*[dx,dy];  % (en mm)
        dup = all(diff(pts,1,1)==0,2); pts(dup,:) = [];
        if size(pts,1) > 2 && ~isequal(pts(1,:),pts(end,:)), pts(end+1,:) = pts(1,:); end
        if size(pts,1) < 4, continue; end
        pts = fftSmooth(pts, K_keep, downsample);
        if polyarea(pts(:,1),pts(:,2)) < minArea_mm2, continue; end
        PS = polyshape(pts(:,1),pts(:,2));
        PS = rmholes(PS);
        PS = simplify(PS, 'KeepCollinear',true);
        % ---- Correction robustesse polyshape
        if PS.NumRegions == 0 || size(PS.Vertices,1) < 4
            warning('Slice %d : polyshape vide', sliceNum); continue;
        end
        [xv, yv] = boundary(PS);
        P = [xv(:), yv(:)];
        dP = diff(P,1,1); dup2 = all(dP==0,2); P(dup2,:) = [];
        if isequal(P(1,:),P(end,:)), P(end,:)=[]; end
        if size(P,1) < 3, continue; end
        cleanContours{end+1} = P;
    end

    if isempty(cleanContours)
        warning('Slice %d : aucun contour valide → skip', sliceNum); continue;
    end

    % === Construction géométrie PDE Toolbox ===
    nS = numel(cleanContours); nv = cellfun(@(P) size(P,1), cleanContours);
    kmax = max(nv); gd = zeros(2+2*kmax, nS); names = cell(1,nS);
    for i = 1:nS
        P = cleanContours{i}; k = size(P,1);
        gd(1,i) = 2; gd(2,i) = k;
        gd(3:2+k, i)   = P(:,1);
        gd(3+k:2+2*k,i)= P(:,2);
        names{i} = sprintf('P%d',i);
    end
    sf = strjoin(names,'+'); ns = char(names)';

    % === Génération mesh PDE ===
    model = createpde();
    try
        [dl, ~] = decsg(gd, sf, ns); geometryFromEdges(model, dl);
    catch ME
        warning('decsg échoué slice %d: geometry problem → skip', sliceNum); continue;
    end
    try
        msh = generateMesh(model, 'Hmax',targetSize, 'Hmin',targetSize, ...
            'Hgrad',1, 'GeometricOrder','linear', ...
            'MesherVersion','preR2013a', 'Jiggle', {'on','minimum'});
    catch ME
        warning('generateMesh échoué slice %d: %s', sliceNum, ME.message); continue;
    end
    g = msh.Nodes'; % [N x 2]  
    H = msh.Elements';  % [M x 3]

    % === Calcul de la conductivité depuis la vraie image CT ===
    idx_ct = sliceNum;   % Attention : slices doivent être alignés !
    ct2D = double(squeeze(vol_ct(:,:,idx_ct)));
    ct2D = (ct2D - min(ct2D(:))) / (max(ct2D(:))-min(ct2D(:))); % Normalise [0,1]
    [Xv, Yv] = meshgrid((0:nx-1)*dx, (0:ny-1)*dy);
    F = griddedInterpolant(Xv', Yv', ct2D', 'linear', 'nearest');
    sig = double(F(g(:,1), g(:,2)));

    % === Electrodes auto (bord) ===
    TR = triangulation(H, g);
    Fb = freeBoundary(TR);
    mids = (g(Fb(:,1),:)+g(Fb(:,2),:))/2;
    cent = mean(mids,1);
    theta = atan2(mids(:,2)-cent(2), mids(:,1)-cent(1));
    [~, idxs] = sort(theta); Fb_s = Fb(idxs,:);
    ds = [0; cumsum(hypot(diff(mids(idxs,1)), diff(mids(idxs,2))))];
    perim = ds(end) + norm(mids(idxs(end),:)-mids(idxs(1),:));
    if norm(mids(idxs(end),:)-mids(idxs(1),:)) > 1e-6
        Fb_s = [Fb_s; Fb_s(1,:)]; ds = [ds; ds(end)+norm(mids(idxs(1),:)-mids(idxs(end),:))];
    end
    elfaces = cell(Ne,1); bounds = linspace(0, perim, Ne+1); curr=1;
    for e = 1:Ne
        inds=[]; while curr<=length(ds) && ds(curr)<bounds(e+1), inds(end+1)=curr; curr=curr+1; end
        if isempty(inds) && e>1, inds=curr-1; elseif isempty(inds), inds=1; end
        segs = unique(Fb_s(inds,:), 'rows');
        if size(segs,2)==1, segs = [segs, circshift(segs,-1)]; segs(end,:)=[]; end
        elfaces{e}=segs;
    end
    for e=1:Ne, if isempty(elfaces{e}), elfaces{e}=Fb_s(e,:); end; end

    % === Sauvegarde .mat ===
    save(fullfile(outDir, sprintf('mesh_slice%03d.mat',sliceNum)), ...
        'g','H','elfaces','sig');
    t_elapsed = toc(t0);
    fprintf('→ mesh_slice%03d.mat sauvegardé (%d noeuds) [%.2fs]\n',sliceNum,size(g,1),t_elapsed);

    % === Preview PNG classique ===
    fig = figure('Visible','off');
    subplot(1,2,1);
    imagesc((0:nx-1)*dx, (0:ny-1)*dy, mask2D); colormap(gray); axis image;
    title(sprintf('Mask bin. slice %d',sliceNum));
    hold on;
    for ci=1:numel(cleanContours), plot(cleanContours{ci}(:,1), cleanContours{ci}(:,2), 'r-', 'LineWidth',1); end
    hold off;
    subplot(1,2,2);
    trisurf(H, g(:,1), g(:,2), zeros(size(g,1),1), sig, ...
        'EdgeColor','none','FaceAlpha',0.8);
    axis equal; view(2); colorbar;
    title('Mesh + Conductivité');
    xlabel('X (mm)'); ylabel('Y (mm)');
    saveas(fig, fullfile(outDir, sprintf('preview_slice%03d.png', sliceNum)));
    close(fig);

    % === Preview PNG "Contours + Mesh" seul ===
    fig2 = figure('Visible','off');
    hold on;
    for ci=1:numel(cleanContours)
        plot(cleanContours{ci}(:,1), cleanContours{ci}(:,2), 'r-', 'LineWidth',2);
    end
    triplot(H, g(:,1), g(:,2), 'k-', 'LineWidth', 0.7);
    hold off;
    axis equal; axis off;
    set(gcf, 'Color', 'w');
    title(sprintf('Contours + Mesh slice %d',sliceNum));
    saveas(fig2, fullfile(outDir, sprintf('contours_mesh_slice%03d.png', sliceNum)));
    close(fig2);
end

fprintf('Tous les meshes et conductivités sont prêts dans %s. [%.2fs]\n', outDir, toc(tic_all));

%% ==== Lissage FFT utilitaire ====
function sm = fftSmooth(pts, K, step)
    d  = hypot(diff(pts(:,1)), diff(pts(:,2)));
    s  = [0; cumsum(d)]/sum(d);
    N  = numel(s);
    su = linspace(0,1,N)';
    Xu = interp1(s,pts(:,1),su,'linear');
    Yu = interp1(s,pts(:,2),su,'linear');
    FX = fft(Xu); FY = fft(Yu);
    K  = min(K,floor((N-1)/2));
    mask = false(N,1); mask([1:K+1, N-K+1:N]) = true;
    FX(~mask)=0; FY(~mask)=0;
    sm = unique([real(ifft(FX)), real(ifft(FY))],'rows','stable');
    if step>1, sm=sm(1:step:end,:); end
    if ~isequal(sm(1,:),sm(end,:)), sm(end+1,:)=sm(1,:); end
end
