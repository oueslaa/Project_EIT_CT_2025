% Appel recommandé dans Plot_Problem_Foward.m :
% plot_slice_ct(img, 'ct_brut.png', info, meta.z_slice, displayMode, plotDir, DPI);
% ou sans export :
% plot_slice_ct(img, [], info, meta.z_slice, displayMode);

function plot_slice_ct(img, out_png, info, z_slice, displayMode, out_dir, dpi)
% PLOT_SLICE_CT  Affiche la slice CT avec la bonne convention d'affichage.
%
% Usage:
%   plot_slice_ct(img, out_png, info, z_slice, displayMode)
%   plot_slice_ct(img, out_png, info, z_slice, displayMode, out_dir, dpi)
%
% Inputs
%   img          : matrice 2D (slice CT déjà extraite)
%   out_png      : nom de fichier PNG (ex: 'ct_brut.png') ou [] pour ne pas exporter
%   info         : struct NIfTI (niftiinfo) pour récupérer l'échelle mm
%   z_slice      : indice de la coupe (affiché dans le titre)
%   displayMode  : 'neurological' ou 'radiological'
%   out_dir      : dossier de sortie (optionnel, par défaut = pwd si out_png non vide)
%   dpi          : résolution export (optionnel, défaut 300)
%
% Effet
%   Ouvre une figure nommée 'CT brut' (NumberTitle off), affiche en niveaux de gris,
%   applique la convention displayMode, et exporte en PNG si out_png est fourni.

    if nargin < 5 || isempty(displayMode), displayMode = 'neurological'; end
    if nargin < 6 || isempty(out_dir), out_dir = pwd; end
    if nargin < 7 || isempty(dpi), dpi = 300; end

    % Détermine si on peut afficher en coordonnées mm
    use_mm = exist('mm_extent_for_slice', 'file') == 2 && ~isempty(info);
    if use_mm
        [H, W] = size(img);
        try
            [XL, YL] = mm_extent_for_slice(info, W, H, z_slice);
        catch
            % fallback en pixels si la fonction/les champs ne conviennent pas
            use_mm = false;
        end
    end

    % --- Figure + axes
    fh = figure('Color','w','Name','CT brut','NumberTitle','off');
    ax = axes('Parent',fh); %#ok<LAXES>

    % --- Affichage image
    if use_mm
        imagesc(ax, XL, YL, img);
        set(ax,'YDir','normal');    % Y vers le haut en mm
    else
        imagesc(ax, img);
    end
    axis(ax,'image'); axis(ax,'off'); colormap(ax, gray);

    % --- Titre + convention d'affichage
    title(ax, sprintf('CT brut - slice %d', z_slice), 'FontWeight','bold');
    if exist('apply_display_convention','file') == 2
        apply_display_convention(ax, displayMode);
    else
        if strcmpi(displayMode,'radiological'), set(ax,'XDir','reverse'); end
    end

    % --- Export éventuel
    if ~isempty(out_png)
        if ~exist(out_dir,'dir'), mkdir(out_dir); end
        exportgraphics(fh, fullfile(out_dir, out_png), 'Resolution', dpi);
    end
end
