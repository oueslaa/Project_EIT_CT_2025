function Plot_Problem_Forward_fn(patient_id, z_slice, varargin)
%PLOT_PROBLEM_FORWARD_FN  Exporte toutes les figures pour (patient_id, z_slice) sans afficher.
% Usage:
%   Plot_Problem_Forward_fn('s0011', 320)                       % enregistre PNG, n'affiche pas
%   Plot_Problem_Forward_fn('s0011', 320, 'ShowFigures', true)  % affiche + enregistre
%
% Options Name-Value :
%   'SavePNG'     : logical, def true
%   'DPI'         : entier, def 300
%   'ShowFigures' : logical, def false

% ---------- Options ----------
p = inputParser;
addParameter(p,'SavePNG',true);
addParameter(p,'DPI',300);
addParameter(p,'ShowFigures',false);
parse(p,varargin{:});
SAVE_PNG  = p.Results.SavePNG;
DPI       = p.Results.DPI;
SHOW_FIGS = p.Results.ShowFigures;

addpath('src_anis', genpath('src_anis'));

% Visibilité des figures
if SHOW_FIGS
    set(groot,'defaultFigureVisible','on');
else
    set(groot,'defaultFigureVisible','off');
end

% ---------- Fichiers ----------
rootOut  = fullfile('Outputs', patient_id, sprintf('slice_%03d', z_slice));
plotDir  = fullfile(rootOut, 'plots_forward');
if ~exist(plotDir,'dir'), mkdir(plotDir); end

meshFile = fullfile(rootOut,'mesh',sprintf('mesh_slice%03d.mat', z_slice));
packFile = fullfile(rootOut,'eit_pack.mat');
vizFile  = fullfile(rootOut,'viz_config.mat');
srcFile  = fullfile(rootOut,'data_sources.mat');
elFile   = fullfile(rootOut,'E_electrodes.mat');

assert(isfile(meshFile), 'Mesh file introuvable: %s', meshFile);
assert(isfile(packFile),'Pack forward introuvable: %s', packFile);
assert(isfile(vizFile), 'viz_config.mat introuvable: %s', vizFile);
assert(isfile(srcFile), 'data_sources.mat introuvable: %s', srcFile);
assert(isfile(elFile),  'E_electrodes.mat introuvable: %s', elFile);

% ---------- Chargements ----------
M   = load(meshFile);   % g, H, triGroup, shapes, domain, contour
S   = load(packFile);   % pack forward
VZ  = load(vizFile);    % viz_config
SRC = load(srcFile);    % meta
EL  = load(elFile);     % E, el_centers

viz  = VZ.viz_config;
meta = SRC.meta;
displayMode = viz.displayMode;

% ====== 1) CT brut =======================================================
try
    [CT, info] = read_nifti3D(meta.ct_file);
    img = squeeze(CT(:,:,meta.z_slice));

    fh = figure('Color','w','Name','CT brut','NumberTitle','off','Visible',onoff(SHOW_FIGS));
    imagesc(img); axis image off; colormap(gray);
    title(sprintf('CT slice z=%d', meta.z_slice));
    apply_display(gca, displayMode);
    savefig_if(fh, fullfile(plotDir,'ct_brut.png'), SAVE_PNG, DPI, SHOW_FIGS);
catch ME
    warning('Impossible d''afficher le CT brut (%s).', ME.message);
end

% ====== 2) Slice segmentee (labelmap + contour) ==========================
try
    if ~exist('info','var')
        [~, info] = read_nifti3D(meta.ct_file);
    end
    outPNG = ternary(SAVE_PNG, fullfile(plotDir,'slice_segmentee.png'), []);
    plot_slice_segmentee(meta.seg_dir, meta.z_slice, M.contour, info, outPNG, displayMode);
catch ME
    warning('plot_slice_segmentee a echoue (%s).', ME.message);
end

% ====== 3) Contours AVANT mesh ===========================================
try
    outFile = ternary(SAVE_PNG, fullfile(plotDir,'contours_avant_mesh.png'), []);
    plot_contours_avant_mesh(M.shapes, outFile, M.domain, M.contour, displayMode);
catch ME
    warning('plot_contours_avant_mesh a echoue (%s).', ME.message);
end

% ====== 4) Mesh seul =====================================================
try
    fh = figure('Color','w','Name','Mesh seul','NumberTitle','off','Visible',onoff(SHOW_FIGS));
    triplot(M.H, M.g(:,1), M.g(:,2), 'Color',[0.6 0.6 0.6]); axis equal tight;
    set(gca,'XTick',[],'YTick',[]); title('Mesh (triangles)');
    apply_display(gca, displayMode);
    savefig_if(fh, fullfile(plotDir,'mesh_only.png'), SAVE_PNG, DPI, SHOW_FIGS);
catch
end

% ====== 5) Mesh + electrodes (numérotées) ================================
try
    fh = figure('Color','w','Name','Mesh + electrodes','NumberTitle','off','Visible',onoff(SHOW_FIGS)); hold on;
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
    axis equal tight; set(gca,'XTick',[],'YTick',[]); title('Mesh + electrodes numerotees');
    apply_display(gca, displayMode);
    savefig_if(fh, fullfile(plotDir,'mesh_electrodes_numbered.png'), SAVE_PNG, DPI, SHOW_FIGS);
catch ME
    warning('Mesh + electrodes : %s', ME.message);
end

% ====== 6) Sigma (forward) ===============================================
try
    Cfaces = ensure_face_cdata(S);
    fh = figure('Color','w','Name','sigma (forward)','NumberTitle','off','Visible',onoff(SHOW_FIGS)); hold on;
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
    if isfield(viz,'sigma_clim') && numel(viz.sigma_clim)==2
        caxis(viz.sigma_clim);
    end
    cb=colorbar; ylabel(cb,'\sigma (S/m)');
    title('mesh + electrodes + \sigma forward'); set(gca,'XTick',[],'YTick',[]);
    apply_display(gca, displayMode);
    savefig_if(fh, fullfile(plotDir,'mesh_electrodes_sigma_forward.png'), SAVE_PNG, DPI, SHOW_FIGS);
catch ME
    warning('Plot sigma(forward) impossible (%s).', ME.message);
end

% ====== 7) Tensions aux electrodes (heatmap + profil) ====================
try
    [Vmat, Ne, K] = ensure_Vmat(S);

    vmax = max(abs(Vmat(:)));
    if vmax < 1e-2, scale = 1e3; unit = 'mV'; else, scale = 1; unit = 'V'; end

    fh = figure('Color','w','Name','Tensions aux electrodes (heatmap)','NumberTitle','off','Visible',onoff(SHOW_FIGS));
    imagesc(scale*Vmat); axis tight; set(gca,'YDir','normal');
    xlabel('Motif d''injection (k)'); ylabel('Electrode (i)');
    cb = colorbar; ylabel(cb, ['Tension (' unit ')']);
    title(sprintf('Tensions aux electrodes (Ne=%d, K=%d)', Ne, K));
    savefig_if(fh, fullfile(plotDir,'electrode_voltages_heatmap.png'), SAVE_PNG, DPI, SHOW_FIGS);

    Vcat = Vmat(:);
    x = 1:numel(Vcat);
    fh2 = figure('Color','w','Name','Profil tensions','NumberTitle','off','Visible',onoff(SHOW_FIGS));
    plot(x, scale*Vcat, '-'); grid on; xlim([1 numel(Vcat)]);
    xlabel(sprintf('Electrodes (1..%d)  [Ne=%d, K=%d]', Ne*K, Ne, K));
    ylabel(['Tension (' unit ')']);
    title('Profil des tensions aux electrodes');
    hold on;
    if K > 1
        xs = (1:K-1)*Ne + 0.5;
        for s = xs
            xline(s, ':', 'HandleVisibility','off');
        end
    end
    if K <= 32, xticks(0:Ne:Ne*K); end
    savefig_if(fh2, fullfile(plotDir,'electrode_voltages_profile.png'), SAVE_PNG, DPI, SHOW_FIGS);
catch ME
    warning('Plot des tensions aux electrodes impossible (%s).', ME.message);
end

% ====== 8) Montage 4x4 potentiels u(x) ===================================
try
    if isfield(S,'U_cell') && ~isempty(S.U_cell) && iscell(S.U_cell)
        Uc_all = S.U_cell(:);
        good = cellfun(@(u) ~isempty(u) && all(isfinite(u)), Uc_all);
        Uc = Uc_all(good);
        if ~isempty(Uc)
            Kall = numel(Uc); kk = min(Kall,16);
            allvals = cell2mat(cellfun(@(u) u(:), Uc(1:kk), 'UniformOutput', false));
            lo = prctile(allvals, 1); hi = prctile(allvals, 99);
            if ~(isfinite(lo) && isfinite(hi) && lo < hi)
                lo = min(allvals); hi = max(allvals);
            end
            fh = figure('Color','w','Name','Potentiels u(x) – premiers motifs','NumberTitle','off','Visible',onoff(SHOW_FIGS));
            tiledlayout(4,4,'Padding','compact','TileSpacing','compact');
            for k=1:kk
                nexttile;
                trisurf(S.H, S.g(:,1), S.g(:,2), Uc{k}, 'EdgeColor','none');
                view(2); axis image off; caxis([lo hi]); title(sprintf('Motif %d',k));
            end
            cb=colorbar('Location','eastoutside'); ylabel(cb,'Potentiel (V)');
            savefig_if(fh, fullfile(plotDir,'potentials_montage_4x4.png'), SAVE_PNG, DPI, SHOW_FIGS);
        else
            warning('U_cell est vide ou invalide : montage potentiels saute.');
        end
    end
catch ME
    warning('Montage potentiels impossible (%s).', ME.message);
end

% ====== 9) |E| moyen =====================================================
try
    if isfield(S,'E_mag_tri') && ~isempty(S.E_mag_tri)
        fh = figure('Color','w','Name','|E| moyen (triangles)','NumberTitle','off','Visible',onoff(SHOW_FIGS));
        patch('Faces',S.H,'Vertices',S.g,'FaceVertexCData',S.E_mag_tri, ...
              'FaceColor','flat','EdgeColor','none'); axis equal tight; colorbar;
        title('|E| moyen sur motifs (V/m)'); set(gca,'XTick',[],'YTick',[]);
        apply_display(gca, displayMode);
        savefig_if(fh, fullfile(plotDir,'Efield_mag_tri.png'), SAVE_PNG, DPI, SHOW_FIGS);
    end
catch ME
    warning('Carte |E| impossible (%s).', ME.message);
end

% ====== 10) Champ E (quiver) =============================================
try
    if isfield(S,'gradUx_tri') && ~isempty(S.gradUx_tri) && isfield(S,'gradUy_tri') && ~isempty(S.gradUy_tri)
        Gx = mean(S.g(S.H,1),2); Gy = mean(S.g(S.H,2),2);
        Ex = -S.gradUx_tri; Ey = -S.gradUy_tri; % E = -∇U
        step = max(1, round(numel(Ex)/800)); idx = 1:step:numel(Ex);
        fh = figure('Color','w','Name','Champ E (quiver)','NumberTitle','off','Visible',onoff(SHOW_FIGS)); hold on;
        patch('Faces',S.H,'Vertices',S.g,'FaceColor','none','EdgeColor',[0.85 0.85 0.85]);
        quiver(Gx(idx), Gy(idx), Ex(idx), Ey(idx), 'AutoScale','on');
        axis equal tight; title('Champ electrique moyen (quiver)'); set(gca,'XTick',[],'YTick',[]);
        apply_display(gca, displayMode);
        savefig_if(fh, fullfile(plotDir,'Efield_quiver.png'), SAVE_PNG, DPI, SHOW_FIGS);
    end
catch ME
    warning('Quiver E impossible (%s).', ME.message);
end

% ====== 11) Motifs de courant I ==========================================
try
    I = ensure_I(S);
    Ne = size(I,1); K  = size(I,2);
    sumI = sum(I, 1);
    fprintf('Verif somme(I) par motif (max abs): %.3e A\n', max(abs(sumI)));

    fh = figure('Color','w','Name','Courants injectes — heatmap','NumberTitle','off','Visible',onoff(SHOW_FIGS));
    imagesc(I); axis tight; set(gca,'YDir','normal');
    xlabel('Motif d''injection (k)'); ylabel('Electrode (i)');
    cb = colorbar; ylabel(cb, 'Courant (A)');
    title(sprintf('Motifs de courant injecte (Ne=%d, K=%d, I_{amp}=%.3g A)', Ne, K, S.params.I_amp));
    savefig_if(fh, fullfile(plotDir, 'currents_heatmap.png'), SAVE_PNG, DPI, SHOW_FIGS);

    Icat = I(:);
    x = 1:numel(Icat);
    fh2 = figure('Color','w','Name','Courants injectes','NumberTitle','off','Visible',onoff(SHOW_FIGS));
    plot(x, Icat, '-'); grid on; xlim([1 numel(Icat)]);
    xlabel(sprintf('Electrodes (1..%d)  [Ne=%d, K=%d]', Ne*K, Ne, K));
    ylabel('Courant (A)');
    title('Courants injectes');
    hold on;
    if K > 1
        xs = (1:K-1)*Ne + 0.5;
        for s = xs
            xline(s, ':', 'HandleVisibility','off');
        end
    end
    if K <= 32, xticks(0:Ne:Ne*K); end
    savefig_if(fh2, fullfile(plotDir, 'currents_profile.png'), SAVE_PNG, DPI, SHOW_FIGS);
catch ME
    warning('Plot motifs de courant impossible (%s).', ME.message);
end

% ====== 12) Reciprocite ==================================================
try
    [Vmat, ~, ~] = ensure_Vmat(S);
    I = ensure_I(S);
    Kc = min(size(Vmat,2), size(I,2));
    Vc = Vmat(:,1:Kc); Ic = I(:,1:Kc);
    Smat = Ic.' * Vc;
    rec_err = norm(Smat - Smat.', 'fro') / max(1e-12, norm(Smat,'fro'));
    fprintf('Reciprocity error (||S - S^T||/||S||): %.3e\n', rec_err);

    fh = figure('Color','w','Name','|I^T V - (I^T V)^T|','NumberTitle','off','Visible',onoff(SHOW_FIGS));
    imagesc(abs(Smat - Smat.')); axis image; set(gca,'YDir','normal');
    title('|I^T V - (I^T V)^T|'); colorbar; xlabel('Motif'); ylabel('Motif');
    savefig_if(fh, fullfile(plotDir,'reciprocity_heatmap.png'), SAVE_PNG, DPI, SHOW_FIGS);
catch ME
    warning('Test de reciprocite impossible (%s).', ME.message);
end

% ---------- Helpers locaux ----------
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

    function [Vmat, Ne, K] = ensure_Vmat(S0)
        if isfield(S0,'Vmat') && ~isempty(S0.Vmat)
            Vmat = S0.Vmat;
            if size(Vmat,1) ~= S0.params.Ne, Vmat = Vmat.'; end
        else
            Vraw = [];
            if isfield(S0,'Imeas') && ~isempty(S0.Imeas), Vraw = S0.Imeas;
            elseif isfield(S0,'Imeas_clean') && ~isempty(S0.Imeas_clean), Vraw = S0.Imeas_clean;
            else, error('Aucune mesure disponible pour construire Vmat.');
            end
            Ne = S0.params.Ne; K = floor(numel(Vraw)/Ne);
            Vmat = reshape(Vraw(1:Ne*K), Ne, K);
            return;
        end
        Ne = size(Vmat,1); K = size(Vmat,2);
    end

    function I = ensure_I(S0)
        if isfield(S0,'Ipat') && ~isempty(S0.Ipat)
            I = S0.Ipat;
        else
            I = EITSim.buildTrigPattern(S0.params.Ne, S0.params.I_amp);
        end
        if size(I,1) ~= S0.params.Ne, I = I.'; end
    end

    function Cfaces = ensure_face_cdata(S0)
        Mt = size(S0.H,1);
        C = S0.sigma_tri;
        if numel(C) == Mt
            Cfaces = C(:);
        elseif numel(C) == size(S0.g,1)
            C = C(:);
            Cfaces = (C(S0.H(:,1)) + C(S0.H(:,2)) + C(S0.H(:,3))) / 3;
        else
            error('Taille de sigma_tri (%d) incompatible avec #faces (%d) et #noeuds (%d).', ...
                  numel(C), Mt, size(S0.g,1));
        end
    end

    function vis = onoff(tf)
        if tf, vis = 'on'; else, vis = 'off'; end
    end

    function savefig_if(fh, outpng, doSave, dpi, doShow)
        if doSave
            exportgraphics(fh, outpng, 'Resolution', dpi);
        end
        if ~doShow
            close(fh);
        end
    end
end
