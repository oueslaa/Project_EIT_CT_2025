function visualise_cube(file)

% Load data
data = load_nii(file);
origin = abs(data.hdr.hist.originator(1:3));
img = data.img;
srow_x = data.hdr.hist.srow_x;
srow_y = data.hdr.hist.srow_y;
srow_z = data.hdr.hist.srow_z;
A = [srow_x;srow_y;srow_z];
R = A(1:3,1:3);
% Obtain voxel coordinates
[Nx,Ny,Nz] = size(img);
i = 1:Nx; j = 1:Ny; k = 1:Nz;

% Obtain RAS coordinates
[i,j,k] = meshgrid(i,j,k);
X = srow_x(1)*i + srow_x(2)*j + srow_x(3)*k + srow_x(4);
Y = srow_y(1)*i + srow_y(2)*j + srow_y(3)*k + srow_y(4);
Z = srow_z(1)*i + srow_z(2)*j + srow_z(3)*k + srow_z(4);
img = img(:); X = X(:); Y = Y(:);Z = Z(:);


% Filter out non-zero voxels
mask =img>0;
img = img(mask);X  =X(mask); Y = Y(mask); Z = Z(mask);
i = i(mask); j = j(mask); k = k(mask);
n = sum(mask);
figure;
hold on;
for ii = 1:n
    % plot_cube(X(i),Y(i),Z(i),1,img(i))
    plot_cube(i(ii),j(ii),k(ii),A,img(ii))
end
axis equal;
grid on;
view(3)

% function plot_cube(x,y,z,h,c)   
% X = [1 -1 -1 1 1]';
% Y = [1 1 -1 -1 1]';
% Z = [-1 -1 -1 -1 -1]';
% X = x + h/2*[zeros(5,1),X,X,zeros(5,1)];
% Y = y+h/2*[zeros(5,1),Y,Y,zeros(5,1)];
% Z =z +h/2*[-ones(5,1),Z,-Z,ones(5,1)];
% c = double(c) +0*X;
% mesh(X,Y,Z,c,'FaceAlpha',0.5);
% end

function plot_cube(i,j,k,A,c)   
I = [1 -1 -1 1 1]';
J = [1 1 -1 -1 1]';
K = [-1 -1 -1 -1 -1]';
I = i + [zeros(5,1),I,I,zeros(5,1)];
J = j+[zeros(5,1),J,J,zeros(5,1)];
K =k +[-ones(5,1),K,-K,ones(5,1)];
X=reshape((A*[I(:),J(:),K(:),ones(20,1)]')',[size(I),3]);
c = double(c) +0*I;
mesh(X(:,:,1),X(:,:,2),X(:,:,3),c,'FaceAlpha',0.5);
