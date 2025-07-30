function wrap = visualise_alpha(file)
% Load data
data = load_nii(file);
origin = abs(data.hdr.hist.originator(1:3));
img = data.img;
srow_x = data.hdr.hist.srow_x;
srow_y = data.hdr.hist.srow_y;
srow_z = data.hdr.hist.srow_z;
A = [srow_x;srow_y;srow_z];
% Obtain voxel coordinates
[Nx,Ny,Nz] = size(img);
i = 1:Nx; j = 1:Ny; k = 1:Nz;

% Obtain RAS coordinates
[i,j,k] = meshgrid(i,j,k);
X = srow_x(1)*i + srow_x(2)*j + srow_x(3)*k + srow_x(4);
Y = srow_y(1)*i + srow_y(2)*j + srow_y(3)*k + srow_y(4);
Z = srow_z(1)*i + srow_z(2)*j + srow_z(3)*k + srow_z(4);
img = img(:); X = X(:); Y = Y(:);Z = Z(:);

% Use miminum eigenvalue of affine transformation for alpha shape
R = A(1:3,1:3);
a = min(abs(eig(R)));

% Filter out non-zero voxels
mask =img>0;
img = img(mask);X  =X(mask); Y = Y(mask); Z = Z(mask);
n = sum(mask);

if var(double(img)) ~= 0
    warning(sprintf("%s is not a binary image",file) )
end
% Apply alphashape wrapping
wrap = alphaShape([X,Y,Z],a);
figure
plot(wrap)