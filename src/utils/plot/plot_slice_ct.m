function plot_slice_ct(img, out_png, info, z_slice)
% NEW signature: passe info et z_slice
[H, W] = size(img);
[XL, YL] = mm_extent_for_slice(info, W, H, z_slice);

fig = figure('Color','w');
imagesc(XL, YL, img);                 % axes en mm
set(gca,'YDir','normal');             % Y vers le haut (comme ton mesh)
axis image; colormap(gray); axis on; box on;
title('CT slice brute','FontSize',18,'FontWeight','bold');

if nargin > 2 && ~isempty(out_png)
    exportgraphics(fig, out_png, 'Resolution', 220);
end
end
