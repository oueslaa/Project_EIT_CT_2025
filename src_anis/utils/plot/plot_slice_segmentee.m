function plot_slice_segmentee(seg_dir, z_slice, ext_smooth, info, outpng, mode)
% PLOT_SLICE_SEGMENTEE  Affiche la labelmap 2D (slice z) + contour extérieur.
%
% But
% ----
% Donner une vue rapide de toutes les segmentations présentes dans seg_dir
% sur la coupe "z_slice", superposées (labelmap), avec le contour extérieur.
%
% Entrées
% -------
% seg_dir   : dossier contenant des .nii(.gz) binaires (masques d'organes)
% z_slice   : indice de coupe (1..D)
% ext_smooth: [N x 2] contour extérieur en mm (optionnel mais recommandé)
% info      : header NIfTI du CT (utile pour passer en coordonnées mm)
% outpng    : chemin du PNG de sortie
% mode      : 'radiological' (défaut) ou 'neurological'
%
% Sortie
% ------
% (figure exportée), pas de variable renvoyée.
%
% Notes
% -----
% - Construit une labelmap en empilant les masques non vides trouvés.
% - La colormap "lines" distingue les labels ; l'ordre suit le tri dir().
%
% Exemple
% -------
% plot_slice_segmentee('Data_set/s0011/segmentations', 310, ext, info, 'slice.png');
%
if nargin < 6 || isempty(mode), mode = 'radiological'; end

S = dir(fullfile(seg_dir, '*.nii*'));

% 1) Construire la labelmap
labelmap = []; seg_names = {}; klabel = 0;
for i = 1:numel(S)
    [V3, ~] = read_nifti3D(fullfile(seg_dir, S(i).name));
    if size(V3,3) < z_slice, continue; end
    mask = squeeze(V3(:,:,z_slice))>0;
    if isempty(labelmap) && any(mask(:)), labelmap = zeros(size(mask)); end
    if any(mask(:))
        klabel = klabel+1; labelmap(mask) = klabel;
        seg_names{klabel} = regexprep(S(i).name,'\.nii(\.gz)?$','');
    end
end
if isempty(labelmap), error('Aucun masque a la slice %d', z_slice); end

[H,W] = size(labelmap);
[XL, YL] = mm_extent_for_slice(info, W, H, z_slice);

fig = figure('Color','w'); hold on;
imagesc(XL, YL, labelmap); set(gca,'YDir','normal'); axis image;
colormap(lines(max(klabel,6)));

% contour exterieur en mm
if ~isempty(ext_smooth)
    plot(ext_smooth(:,1), ext_smooth(:,2), 'k-', 'LineWidth', 3);
end

title('Image segmentee + contour ext', 'FontWeight','bold','FontSize',16);

try, hide_axes_toolbar(gca); end %#ok<TRYNC>
apply_display_convention(gca, mode);

exportgraphics(fig, outpng, 'Resolution', 220);
end
