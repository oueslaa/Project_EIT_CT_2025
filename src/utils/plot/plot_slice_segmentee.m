function plot_slice_segmentee(seg_dir, z_slice, ext_smooth, info, outpng, mode)
% Affiche la labelmap de la slice en coordonnées mm + le contour extérieur.
% mode (optionnel): 'radiological' (défaut) ou 'neurological'

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

% Orientation + toolbar
try, hide_axes_toolbar(gca); end %#ok<TRYNC>
apply_display_convention(gca, mode);

exportgraphics(fig, outpng, 'Resolution', 220);
end
