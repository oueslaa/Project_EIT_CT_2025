function [XL, YL] = mm_extent_for_slice(info, W, H, z)
% Renvoie les bornes [xmin xmax], [ymin ymax] en mm pour la slice z
A = info.Transform.T;                 % voxel -> mm
corn = [ 1  1  z;
         W  1  z;
         1  H  z;
         W  H  z ];
corn = (A * [corn, ones(4,1)]')';     % -> mm
XL = [min(corn(:,1)) max(corn(:,1))];
YL = [min(corn(:,2)) max(corn(:,2))];
end
