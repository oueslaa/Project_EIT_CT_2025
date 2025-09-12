
% ============================
% File: Plot_Problem_Forward.m
% ============================
% VISUALISATION DU PROBLÈME FORWARD EIT
% -------------------------------------------------------------------------
% BUT
%   Afficher TOUTES les figures utiles à partir des fichiers produits par
%   Problem_Forward.m : CT brut, contours/mesh, électrodes, champs σ, |E|,
%   potentiels, métriques (SNR, réciprocité), etc.
%
% ENTRÉES (fichiers attendus dans Outputs/<patient>/slice_xxx/)
%   - mesh/mesh_sliceXXX.mat : g, H, triGroup, shapes, domain, contour
%   - eit_pack.mat           : pack forward (Vmat/Imeas/Ipat, sigma_tri, params…)
%   - viz_config.mat         : displayMode, sigma_clim, cmap
%   - data_sources.mat       : meta (ct_file, seg_dir, z_slice, patient_id)
%   - E_electrodes.mat       : E (cell d'arêtes du bord), el_centers
%
% SORTIES
%   - Dossier plots_forward/ avec PNG haute définition (si SAVE_PNG==true)
%
% REMARQUES / PIÈGES
%   - Ce script AFFICHE les figures (contrairement à Problem_Forward).
%   - displayMode:
%       * 'radiological'  : convention DICOM (X inversé visuellement)
%       * 'neurological'  : X normal
%   - Les figures utilisent des garde-fous (try/catch) pour ne pas interrompre
%     tout le script si un plot échoue.
% -------------------------------------------------------------------------

% (Laisse ce fichier en script; toutes les fonctions locales sont à la fin)

clear; close all; clc;
addpath('src', genpath('src'));

% --- Ici on veut VOIR les figures ----------------------------------------
set(groot,'defaultFigureVisible','on');

% --- Sauvegarde PNG ? ----------------------------------------------------
SAVE_PNG = true;         % true -> export des PNG
DPI      = 300;          % résolution des exports

% --- Paramètres utilisateur ----------------------------------------------
patient_id  = 's0011';
z_slice     = 301;

rootOut     = fullfile('Outputs', patient_id, sprintf('slice_%03d', z_slice));
plotDir     = fullfile(rootOut, 'plots_forward');
if ~exist(plotDir,'dir'), mkdir(plotDir); end

meshFile    = fullfile(rootOut, 'mesh', sprintf('mesh_slice%03d.mat', z_slice));
packFile    = fullfile(rootOut, 'eit_pack.mat');
vizFile     = fullfile(rootOut, 'viz_config.mat');
srcFile     = fullfile(rootOut, 'data_sources.mat');
elFile      = fullfile(rootOut, 'E_electrodes.mat');

% --- Vérifications de présence des fichiers ------------------------------
assert(isfile(meshFile), 'Mesh file introuvable: %s', meshFile);
assert(isfile(packFile), 'Pack forward introuvable: %s', packFile);
assert(isfile(vizFile),  'viz_config.mat introuvable: %s', vizFile);
assert(isfile(srcFile),  'data_sources.mat introuvable: %s', srcFile);
assert(isfile(elFile),   'E_electrodes.mat introuvable: %s', elFile);

% --- Chargements ----------------------------------------------------------
M    = load(meshFile);          % g, H, triGroup, shapes, domain, contour
S    = load(packFile);          % pack forward complet
VZ   = load(vizFile);           % viz_config
SRC  = load(srcFile);           % meta
EL   = load(elFile);            % E, el_centers

viz  = VZ.viz_config;
meta = SRC.meta;
displayMode = viz.displayMode;

%% 1) CT brut (robuste) ----------------------------------------------------
try
    [CT, info] = read_nifti3D(meta.ct_file);
    img = squeeze(CT(:,:,meta.z_slice));

    fh = figure('Color','w','Name','CT brut','NumberTitle','off');
    imagesc(img); axis image off; colormap(gray);
    title(sprintf('CT slice z=%d', meta.z_slice));
    apply_display(gca, displayMode);
    if SAVE_PNG
        exportgraphics(fh, fullfile(plotDir,'ct_brut.png'), 'Resolution', DPI);
    end
catch ME
    warning('Impossible d''afficher le CT brut (%s).', ME.message);
end

%% 2) Slice segmentée (labelmap + contour extérieur) ----------------------
try
    if ~exist('info','var')
        [~, info] = read_nifti3D(meta.ct_file);
    end
    outPNG = ternary(SAVE_PNG, fullfile(plotDir,'slice_segmentee.png'), []);
    plot_slice_segmentee(meta.seg_dir, meta.z_slice, M.contour, info, outPNG, displayMode);
catch ME
    warning('plot_slice_segmentee a echoué (%s).', ME.message);
end

%% 3) Contours AVANT mesh --------------------------------------------------
try
    outFile = ternary(SAVE_PNG, fullfile(plotDir, 'contours_avant_mesh.png'), []);
    plot_contours_avant_mesh(M.shapes, outFile, M.domain, M.contour, displayMode);
catch ME
    warning('plot_contours_avant_mesh a echoué (%s).', ME.message);
end

%% 4) Mesh seul ------------------------------------------------------------
try
    fh = figure('Color','w','Name','Mesh seul','NumberTitle','off');
    triplot(M.H, M.g(:,1), M.g(:,2), 'Color',[0.6 0.6 0.6]); axis equal tight;
    set(gca,'XTick',[],'YTick',[]); title('Mesh (triangles)');
    apply_display(gca, displayMode);
    if SAVE_PNG, exportgraphics(fh, fullfile(plotDir,'mesh_only.png'), 'Resolution', DPI); end
catch
    % rien
end

%% 5) Mesh + électrodes (numérotées) --------------------------------------
try
    fh = figure('Color','w','Name','Mesh + électrodes','NumberTitle','off'); hold on;
    triplot(M.H, M.g(:,1), M.g(:,2), 'Color',[0.8 0.8 0.8]);
    for k = 1:S.params.Ne
        Ek = EL.E{k};
        plot(M.g(Ek(:,1),1), M.g(Ek(:,1),2), '-', 'Color', [0.2 0.2 1], 'LineWidth', 2);
    end
    if ~isempty(EL.el_centers)
        scatter(EL.el_centers(:,1), EL.el_centers(:,2), 50, 'b', 'filled');
        for k = 1:size(EL.el_centers,1)
            text(EL.el_centers(k,1), EL.el_centers(k,2), sprintf('%d',k), ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'Color','w','FontWeight','bold');
        end
    end
    axis equal tight; set(gca,'XTick',[],'YTick',[]); title('Mesh + électrodes numérotées');
    apply_display(gca, displayMode);
    if SAVE_PNG, exportgraphics(fh, fullfile(plotDir,'mesh_electrodes_numbered.png'), 'Resolution', DPI); end
catch ME
    warning('Mesh + électrodes : %s', ME.message);
end

%% 6) Sigma (forward) — mesh + électrodes ---------------------------------
try
    Cfaces = ensure_face_cdata(S);   % garantit 1 valeur par triangle
    fh = figure('Color','w','Name','sigma (forward)','NumberTitle','off'); hold on;
    patch('Faces',S.H,'Vertices',S.g,'FaceVertexCData',Cfaces, ...
          'FaceColor','flat','EdgeColor',[0.75 0.75 0.75]);
    for k = 1:S.params.Ne
        Ek = EL.E{k};
        plot(S.g(Ek(:,1),1), S.g(Ek(:,1),2), '-', 'Color', [0.2 0.2 1], 'LineWidth', 2.5);
    end
    if ~isempty(EL.el_centers)
        scatter(EL.el_centers(:,1), EL.el_centers(:,2), 60, 'b', 'filled');
    end
    axis equal tight; colormap(turbo);
    if isfield(viz,'sigma_clim') && all(isfinite(viz.sigma_clim)) && numel(viz.sigma_clim)==2
        caxis(viz.sigma_clim);
    end
    cb=colorbar; ylabel(cb,'\sigma (S/m)');
    title('mesh + electrodes + \sigma forward'); set(gca,'XTick',[],'YTick',[]);
    apply_display(gca, displayMode);
    if SAVE_PNG
        exportgraphics(fh, fullfile(plotDir,'mesh_electrodes_sigma_forward.png'), 'Resolution', DPI);
    end
catch ME
    warning('Plot sigma(forward) impossible (%s).', ME.message);
end

%% 7) Tensions aux électrodes (heatmap + profil ) -------------------------
try
    [Vmat, Ne, K] = ensure_Vmat(S);

    % --- unité lisible ---
    vmax = max(abs(Vmat(:)));
    if vmax < 1e-2, scale = 1e3; unit = 'mV'; else, scale = 1; unit = 'V'; end

    % --- Heatmap (Ne x K) ---
    fh = figure('Color','w','Name','Tensions aux electrodes (heatmap)','NumberTitle','off');
    imagesc(scale*Vmat); axis tight; set(gca,'YDir','normal');
    xlabel('Motif d''injection (k)'); ylabel('Electrode (i)');
    cb = colorbar; ylabel(cb, ['Tension (' unit ')']);
    title(sprintf('Tensions aux electrodes (Ne=%d, K=%d)', Ne, K));
    if SAVE_PNG
        exportgraphics(fh, fullfile(plotDir,'electrode_voltages_heatmap.png'), 'Resolution', DPI);
    end

    % --- Profil sur 1..Ne*K ---
    Vcat = Vmat(:);
    x = 1:numel(Vcat);

    fh2 = figure('Color','w','Name','Profil tensions','NumberTitle','off');
    plot(x, scale*Vcat, '-'); grid on; xlim([1 numel(Vcat)]);
    xlabel(sprintf('Électrodes (1..%d)  [Ne=%d, K=%d]', Ne*K, Ne, K));
    ylabel(sprintf('Tension (%s)', unit));
    title('Profil des tensions aux électrodes');

    % séparateurs à chaque frontière de motif
    hold on;
    if K > 1
        xs = (1:K-1)*Ne + 0.5;
        for s = xs
            xline(s, ':', 'HandleVisibility','off');
        end
    end
    if K <= 32, xticks(0:Ne:Ne*K); end

    if SAVE_PNG
        exportgraphics(fh2, fullfile(plotDir,'electrode_voltages_profile.png'), 'Resolution', DPI);
    end
catch ME
    warning('Plot des tensions aux electrodes impossible (%s).', ME.message);
end

%% 8) Montage 4x4 des potentiels nodaux u(x) -------------------------------
try
    if isfield(S,'U_cell') && ~isempty(S.U_cell) && iscell(S.U_cell)
        Uc_all = S.U_cell(:);
        good = cellfun(@(u) ~isempty(u) && all(isfinite(u)), Uc_all);
        Uc = Uc_all(good);
        if ~isempty(Uc)
            Kall = numel(Uc); kk = min(Kall,16);

            % caxis global robuste via percentiles (évite outliers)
            allvals = cell2mat(cellfun(@(u) u(:), Uc(1:kk), 'UniformOutput', false));
            lo = prctile(allvals, 1); hi = prctile(allvals, 99);
            if ~(isfinite(lo) && isfinite(hi) && lo < hi), lo = min(allvals); hi = max(allvals); end

            fh = figure('Color','w','Name','Potentiels u(x) – premiers motifs','NumberTitle','off');
            tiledlayout(4,4,'Padding','compact','TileSpacing','compact');
            for k=1:kk
                nexttile;
                trisurf(S.H, S.g(:,1), S.g(:,2), Uc{k}, 'EdgeColor','none');
                view(2); axis image off; caxis([lo hi]); title(sprintf('Motif %d',k));
            end
            cb=colorbar('Location','eastoutside'); ylabel(cb,'Potentiel (V)');
            if SAVE_PNG
                exportgraphics(fh, fullfile(plotDir,'potentials_montage_4x4.png'), 'Resolution', DPI);
            end
        else
            warning('U_cell est vide ou invalide : montage potentiels sauté.');
        end
    end
catch ME
    warning('Montage potentiels impossible (%s).', ME.message);
end

%% 9) Champ électrique |E| moyen sur motifs --------------------------------
try
    if isfield(S,'E_mag_tri') && ~isempty(S.E_mag_tri)
        fh = figure('Color','w','Name','|E| moyen (triangles)','NumberTitle','off');
        patch('Faces',S.H,'Vertices',S.g,'FaceVertexCData',S.E_mag_tri, ...
              'FaceColor','flat','EdgeColor','none'); axis equal tight; colorbar;
        title('|E| moyen sur motifs (V/m)'); set(gca,'XTick',[],'YTick',[]);
        apply_display(gca, displayMode);
        if SAVE_PNG, exportgraphics(fh, fullfile(plotDir,'Efield_mag_tri.png'), 'Resolution', DPI); end
    end
catch ME
    warning('Carte |E| impossible (%s).', ME.message);
end

%% 10) Champ E (quiver sur barycentres) -----------------------------------
try
    if isfield(S,'gradUx_tri') && ~isempty(S.gradUx_tri) && isfield(S,'gradUy_tri') && ~isempty(S.gradUy_tri)
        % Barycentres de triangles
        Gx = mean(S.g(S.H,1),2); Gy = mean(S.g(S.H,2),2);
        Ex = -S.gradUx_tri; Ey = -S.gradUy_tri; % E = -∇U

        % Sous-échantillonnage pour lisibilité
        step = max(1, round(numel(Ex)/800)); idx = 1:step:numel(Ex);

        fh = figure('Color','w','Name','Champ E (quiver)','NumberTitle','off'); hold on;
        patch('Faces',S.H,'Vertices',S.g,'FaceColor','none','EdgeColor',[0.85 0.85 0.85]);
        quiver(Gx(idx), Gy(idx), Ex(idx), Ey(idx), 'AutoScale','on');
        axis equal tight; title('Champ électrique moyen (quiver)'); set(gca,'XTick',[],'YTick',[]);
        apply_display(gca, displayMode);
        if SAVE_PNG, exportgraphics(fh, fullfile(plotDir,'Efield_quiver.png'), 'Resolution', DPI); end
    end
catch ME
    warning('Quiver E impossible (%s).', ME.message);
end

%% 11) Motifs de courant injecté (I) — heatmap + profils ------------------
try
    I = ensure_I(S);
    Ne = size(I,1); K  = size(I,2);

    % --- vérif somme nulle par motif ---
    sumI = sum(I, 1);
    fprintf('Vérif somme(I) par motif (max abs): %.3e A\n', max(abs(sumI)));

    % --- Heatmap ---
    fh = figure('Color','w','Name','Courants injectés — heatmap','NumberTitle','off');
    imagesc(I); axis tight; set(gca,'YDir','normal');
    xlabel('Motif d''injection (k)'); ylabel('Électrode (i)');
    cb = colorbar; ylabel(cb, 'Courant (A)');
    title(sprintf('Motifs de courant injecté (Ne=%d, K=%d, I_{amp}=%.3g A)', Ne, K, S.params.I_amp));
    if SAVE_PNG
        exportgraphics(fh, fullfile(plotDir, 'currents_heatmap.png'), 'Resolution', DPI);
    end

    % --- Profil sur 1..Ne*K ---
    Icat = I(:);
    x = 1:numel(Icat);

    fh2 = figure('Color','w','Name','Courants injectés','NumberTitle','off');
    plot(x, Icat, '-'); grid on; xlim([1 numel(Icat)]);
    xlabel(sprintf('Électrodes (1..%d)  [Ne=%d, K=%d]', Ne*K, Ne, K));
    ylabel('Courant (A)');
    title('Courants injectés');

    % Séparateurs entre motifs
    hold on;
    if K > 1
        xs = (1:K-1)*Ne + 0.5;  % entre deux blocs d'électrodes
        for s = xs
            xline(s, ':', 'HandleVisibility','off');
        end
    end
    if K <= 32, xticks(0:Ne:Ne*K); end

    if SAVE_PNG
        exportgraphics(fh2, fullfile(plotDir, 'currents_profile.png'), 'Resolution', DPI);
    end
catch ME
    warning('Plot motifs de courant impossible (%s).', ME.message);
end

%% 12) Réciprocité (||I^T V - (I^T V)^T||) --------------------------------
try
    [Vmat, ~, ~] = ensure_Vmat(S);
    I = ensure_I(S);

    % Aligne sur le plus petit K commun
    Kc = min(size(Vmat,2), size(I,2));
    Vc = Vmat(:,1:Kc); Ic = I(:,1:Kc);

    Smat = Ic.' * Vc;  % K x K (symétrique si réciprocité parfaite)
    rec_err = norm(Smat - Smat.', 'fro') / max(1e-12, norm(Smat,'fro'));
    fprintf('Reciprocity error (||S - S^T||/||S||): %.3e\n', rec_err);

    fh = figure('Color','w','Name','|I^T V - (I^T V)^T|','NumberTitle','off');
    imagesc(abs(Smat - Smat.')); axis image; set(gca,'YDir','normal');
    title('|I^T V - (I^T V)^T|'); colorbar; xlabel('Motif'); ylabel('Motif');
    if SAVE_PNG
        exportgraphics(fh, fullfile(plotDir,'reciprocity_heatmap.png'), 'Resolution', DPI);
    end
catch ME
    warning('Test de réciprocité impossible (%s).', ME.message);
end

%% 13) Qualité de mesh: aires et angle minimal -----------------------------
try
    P1 = S.g(S.H(:,1),:); P2 = S.g(S.H(:,2),:); P3 = S.g(S.H(:,3),:);
    A  = 0.5*abs( P1(:,1).*(P2(:,2)-P3(:,2)) + P2(:,1).*(P3(:,2)-P1(:,2)) + P3(:,1).*(P1(:,2)-P2(:,2)) );
    fh = figure('Color','w','Name','Histogramme aires triangles (mm^2)','NumberTitle','off');
    histogram(A, 40); xlabel('Aire (mm^2)'); ylabel('Comptes'); grid on;
    if SAVE_PNG, exportgraphics(fh, fullfile(plotDir, 'mesh_area_hist.png'), 'Resolution', DPI); end

    v1 = P2-P1; v2 = P3-P2; v3 = P1-P3;
    ang = @(u,v) acosd( max(-1,min(1, sum(u.*v,2)./(sqrt(sum(u.^2,2)).*sqrt(sum(v.^2,2))))) );
    a1 = ang(v1,-v3); a2 = ang(v2,-v1); a3 = ang(v3,-v2);
    amin = min([a1 a2 a3],[],2);
    fh = figure('Color','w','Name','Histogramme angle minimal (deg)','NumberTitle','off');
    histogram(amin, 40); xlabel('Angle minimal (°)'); ylabel('Comptes'); grid on;
    if SAVE_PNG, exportgraphics(fh, fullfile(plotDir, 'mesh_minangle_hist.png'), 'Resolution', DPI); end
catch ME
    warning('Qualité de mesh non tracée (%s).', ME.message);
end

%% 14) Histogramme des conductivités élémentaires --------------------------
try
    fh = figure('Color','w','Name','Histogramme \\sigma (par triangle)','NumberTitle','off');
    histogram(S.sigma_tri, 30); xlabel('\\sigma (S/m)'); ylabel('Comptes'); grid on;
    if isfield(viz,'sigma_clim') && numel(viz.sigma_clim)==2
        xlim(viz.sigma_clim);
    end
    if SAVE_PNG, exportgraphics(fh, fullfile(plotDir, 'sigma_hist.png'), 'Resolution', DPI); end
catch ME
    warning('Histogramme sigma impossible (%s).', ME.message);
end

%% 15) SNR et distribution du bruit ---------------------------------------
try
    assert(isfield(S,'Imeas_clean') && ~isempty(S.Imeas_clean) && isfield(S,'Imeas') && ~isempty(S.Imeas), ...
        'Imeas_clean/Imeas manquants pour SNR.');
    noise = S.Imeas - S.Imeas_clean;
    snr = 20*log10( norm(S.Imeas_clean) / max(1e-12, norm(noise)) );
    fprintf('SNR (dB) ~ %.2f dB\n', snr);
    fh = figure('Color','w','Name','Histogramme bruit','NumberTitle','off');
    histogram(noise, 40); xlabel('Bruit (V)'); ylabel('Comptes'); grid on;
    if SAVE_PNG, exportgraphics(fh, fullfile(plotDir, 'noise_hist.png'), 'Resolution', DPI); end
catch ME
    warning('SNR/bruit non tracés (%s).', ME.message);
end

% ===================== Helpers locaux (FONCTIONS EN FIN DE SCRIPT) =====================

function out = ternary(cond,a,b)
    if cond, out = a; else, out = b; end
end

function apply_display(ax, mode)
    if strcmpi(mode,'radiological')
        set(ax,'XDir','reverse');
    else
        set(ax,'XDir','normal');
    end
end

function [Vmat, Ne, K] = ensure_Vmat(S)
    % Garantit Vmat de taille Ne x K à partir de S (pack forward)
    if isfield(S,'Vmat') && ~isempty(S.Vmat)
        Vmat = S.Vmat;
        if size(Vmat,1) ~= S.params.Ne, Vmat = Vmat.'; end
    else
        Vraw = [];
        if isfield(S,'Imeas') && ~isempty(S.Imeas), Vraw = S.Imeas;
        elseif isfield(S,'Imeas_clean') && ~isempty(S.Imeas_clean), Vraw = S.Imeas_clean;
        else, error('Aucune mesure disponible pour construire Vmat.');
        end
        Ne = S.params.Ne; K = floor(numel(Vraw)/Ne);
        Vmat = reshape(Vraw(1:Ne*K), Ne, K);
        return;
    end
    Ne = size(Vmat,1); K = size(Vmat,2);
end

function I = ensure_I(S)
    if isfield(S,'Ipat') && ~isempty(S.Ipat)
        I = S.Ipat;
    else
        I = EITSim.buildTrigPattern(S.params.Ne, S.params.I_amp);
    end
    if size(I,1) ~= S.params.Ne, I = I.'; end
end

function Cfaces = ensure_face_cdata(S)
    Mt = size(S.H,1);
    C = S.sigma_tri;
    if numel(C) == Mt
        Cfaces = C(:);
    elseif numel(C) == size(S.g,1)
        C = C(:);
        Cfaces = (C(S.H(:,1)) + C(S.H(:,2)) + C(S.H(:,3))) / 3;
    else
        error('Taille de sigma_tri (%d) incompatible avec #faces (%d) et #noeuds (%d).',...
              numel(C), Mt, size(S.g,1));
    end
end
