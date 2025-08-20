% ============================
% File: Plot_Problem_Foward.m
% ============================
% Affiche TOUT le problème direct EIT à partir des fichiers produits par
% Problem_Foward (qui reste batch-only). Ce script:
%  - affiche CT brut
%  - affiche la slice segmentee (fonction d'origine)
%  - affiche les "contours avant mesh"
%  - affiche le mesh seul / mesh + electrodes
%  - affiche sigma (forward) sur mesh + electrodes
%  - (optionnel) appelle simulate_inverse([],true) pour afficher les plots internes

clear; close all; clc;
addpath('src', genpath('src'));

% --- Ici on veut VOIR les figures
set(groot,'defaultFigureVisible','on');

% --- Sauvegarde PNG ?
SAVE_PNG = true;
DPI      = 300;

% --- Paramètres utilisateur ---
patient_id  = 's0011';
z_slice     = 202;

rootOut     = fullfile('Outputs', patient_id, sprintf('slice_%03d', z_slice));
plotDir     = fullfile(rootOut, 'plots_forward');
if ~exist(plotDir,'dir'), mkdir(plotDir); end

meshFile    = fullfile(rootOut, 'mesh', sprintf('mesh_slice%03d.mat', z_slice));
packFile    = fullfile(rootOut, 'eit_pack.mat');
vizFile     = fullfile(rootOut, 'viz_config.mat');
srcFile     = fullfile(rootOut, 'data_sources.mat');
elFile      = fullfile(rootOut, 'E_electrodes.mat');

% --- Vérifs ---
assert(isfile(meshFile), 'Mesh file introuvable: %s', meshFile);
assert(isfile(packFile), 'Pack forward introuvable: %s', packFile);
assert(isfile(vizFile),  'viz_config.mat introuvable: %s', vizFile);
assert(isfile(srcFile),  'data_sources.mat introuvable: %s', srcFile);
assert(isfile(elFile),   'E_electrodes.mat introuvable: %s', elFile);

M    = load(meshFile);          % g, H, triGroup, shapes, domain, contour
S    = load(packFile);          % g,H,triGroup,E,el_centers,sigma_tri,params,Imeas,domain
VZ   = load(vizFile);           % viz_config
SRC  = load(srcFile);           % meta
EL   = load(elFile);            % E, el_centers
viz  = VZ.viz_config;
meta = SRC.meta;

displayMode = viz.displayMode;

%% 1) CT brut

try
    [CT, info] = read_nifti3D(meta.ct_file);
    img = squeeze(CT(:,:,meta.z_slice));
    if SAVE_PNG
        plot_slice_ct(img, 'ct_brut.png', info, meta.z_slice, displayMode, plotDir, DPI);
    else
        plot_slice_ct(img, [], info, meta.z_slice, displayMode);
    end
catch ME
    warning('Impossible d''afficher le CT brut (%s).', ME.message);
end

%% 2) Slice segmentee — fonction d'origine
try
    if ~exist('info','var')
        [~, info] = read_nifti3D(meta.ct_file);
    end
    if SAVE_PNG, outPNG = fullfile(plotDir, 'slice_segmentee.png'); else, outPNG = []; end
    plot_slice_segmentee(meta.seg_dir, meta.z_slice, M.contour, info, outPNG, displayMode);
catch ME
    warning('plot_slice_segmentee a echoue (%s).', ME.message);
end

%% 3) Contours AVANT mesh (organes visibles avant la triangulation)
%   -> on réutilise shapes/domain/contour sauvegardés par MeshBuilder
try
    fh = figure('Color','w','Name','Contours avant mesh','NumberTitle','off');
    % plot_contours_avant_mesh(shapes, outFile, domain, ext_smooth, displayMode)
    if SAVE_PNG
        outFile = fullfile(plotDir, 'contours_avant_mesh.png');
    else
        outFile = [];
    end
    % Beaucoup de versions de cette fonction acceptent outFile=[] pour juste afficher.
    % Si la tienne exige un chemin, on lui passe outFile; sinon, elle ignorera.
    plot_contours_avant_mesh(M.shapes, outFile, M.domain, M.contour, displayMode);
    if ~SAVE_PNG
        % rien à faire
    end
catch ME
    warning('plot_contours_avant_mesh a echoue (%s).', ME.message);
end


%% 6) Sigma (forward) — mesh + electrodes 
try
    fh = figure('Color','w','Name','sigma (forward)','NumberTitle','off'); hold on;
    patch('Faces',S.H,'Vertices',S.g,'FaceVertexCData',S.sigma_tri, ...
          'FaceColor','flat','EdgeColor',[0.7 0.7 0.7]);
    triplot(S.H, S.g(:,1), S.g(:,2), 'Color', [0.6 0.6 0.6]);
    for k = 1:S.params.Ne
        Ek = EL.E{k};
        for e = 1:size(Ek,1)
            plot(S.g(Ek(e,:),1), S.g(Ek(e,:),2), '-', 'Color', [0.2 0.2 1], 'LineWidth', 3);
        end
    end
    if ~isempty(EL.el_centers)
        scatter(EL.el_centers(:,1), EL.el_centers(:,2), 70, 'b', 'filled');
    end
    axis equal tight;
    colormap(turbo); 
    caxis(viz.sigma_clim);
    cb=colorbar; ylabel(cb,'S/m');
    title('mesh + electrodes + sigma forward ','FontWeight','bold'); set(gca,'XTick',[],'YTick',[]);
    apply_display_convention(gca, displayMode);
    if SAVE_PNG
        exportgraphics(fh, fullfile(plotDir, 'mesh_electrodes_sigma_forward.png'), 'Resolution', DPI);
    end
catch ME
    warning('Plot sigma(forward) impossible (%s).', ME.message);
end

%% 7) Tensions aux electrodes (forward)
try
    V = S.Imeas;
    assert(~isempty(V), 'S.Imeas est vide. As-tu bien lance Problem_Foward ?');
    Ne = S.params.Ne;

    % On travaille sur la partie reelle (ou le module si complexe)
    if ~isreal(V), Vp = abs(V); else, Vp = V; end

    % Tentative de remise en forme: [Ne x K] si possible
    n  = numel(Vp);
    kf = n / Ne;
    is_matrix = abs(kf - round(kf)) < 1e-12;
    if is_matrix
        K = round(kf);
        Vmat = reshape(Vp, [Ne, K]);
    else
        Vmat = Vp(:);  % restera un vecteur
    end

    % Choix d'unite: V ou mV (auto)
    vmax = max(abs(Vp));
    if vmax < 1e-2
        scale = 1e3; unit = 'mV';
    else
        scale = 1;   unit = 'V';
    end

    if is_matrix
        % ---- Heatmap electrodes x motifs d'injection ----
        fh = figure('Color','w','Name','Tensions aux electrodes (heatmap)','NumberTitle','off');
        imagesc(1:K, 1:Ne, scale*Vmat);
        axis tight; set(gca,'YDir','normal');
        xlabel('Motif d''injection'); ylabel('Electrode');
        cb = colorbar; ylabel(cb, ['Tension (' unit ')']);
        title('Tensions aux electrodes (forward)');
        if SAVE_PNG
            exportgraphics(fh, fullfile(plotDir, 'electrode_voltages_heatmap.png'), 'Resolution', DPI);
        end

        % ---- Profil pour le premier motif ----
        fh2 = figure('Color','w','Name','Profil tensions - motif #1','NumberTitle','off');
        plot(1:Ne, scale*Vmat(:,1), '-o'); grid on;
        xlim([1 Ne]); xticks(1:Ne);
        xlabel('Electrode'); ylabel(['Tension (' unit ')']);
        title('Profil des tensions — motif #1');
        if SAVE_PNG
            exportgraphics(fh2, fullfile(plotDir, 'electrode_voltages_profile_motif01.png'), 'Resolution', DPI);
        end
    else
        % ---- Simple profil si on ne peut pas reshaper ----
        fh = figure('Color','w','Name','Tensions aux electrodes (profil)','NumberTitle','off');
        plot(1:numel(Vmat), scale*Vmat, '-o'); grid on;
        xlabel('Index mesure'); ylabel(['Tension (' unit ')']);
        title('Tensions aux electrodes (forward)');
        if SAVE_PNG
            exportgraphics(fh, fullfile(plotDir, 'electrode_voltages_profile.png'), 'Resolution', DPI);
        end
    end
catch ME
    warning('Plot des tensions aux electrodes impossible (%s).', ME.message);
end



% 1) somme de courants
inj = EITSim.buildTrigPattern(S.params.Ne, S.params.I_amp);
assert(all(abs(sum(inj,1))<1e-12))

% 2) réciprocité (si Vmat = Ne x K)

V = reshape(S.Imeas_clean, S.params.Ne, []);   % Ne x K
I = EITSim.buildTrigPattern(S.params.Ne, S.params.I_amp);   % Ne x K (mêmes drives)

Smat = I.' * V;                                 % K x K
recip_err = norm(Smat - Smat.', 'fro') / norm(Smat, 'fro');
fprintf('Reciprocity error (symmetry of I^T V): %.3e\n', recip_err);



% 3) corrélation colonne k et k+1 décalée :
corrs = zeros(1,K-1);
for k=1:K-1, corrs(k) = max(xcorr(Vmat(:,k),circshift(Vmat(:,k+1),1),'coeff')); end
fprintf('Corrélation moyenne décalage: %.3f\n', mean(corrs))





%% --------- Utilitaires ----------
function apply_display_convention(ax, mode)
if nargin<1 || isempty(ax), ax=gca; end
if strcmpi(mode,'radiological')
    set(ax,'XDir','reverse');
else
    set(ax,'XDir','normal');
end
end
