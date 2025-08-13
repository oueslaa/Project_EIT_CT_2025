%% main_pipeline.m - Pipeline EIT 2D slice par slice (refactor)
clear; close all; clc;
addpath('src', genpath('src'));
addpath(genpath('/Users/anis/Documents/StageInria/Code/OOEIT-main'));

% ==== Figures visibles ====
set(groot,'defaultFigureVisible','on');

% ==== Paramètres utilisateur ====
patient_id = 's0011';
z_slice    = 301;
rootOut    = fullfile('Outputs', patient_id, sprintf('slice_%03d', z_slice));
if ~exist(rootOut,'dir'), mkdir(rootOut); end

% Convention d'affichage (cohérente partout)
displayMode = 'neurological';  % ou 'radiological'
PLOT = struct('dpi',300,'pdf',true,'svg',false,'fig',false,'mode',displayMode); %#ok<NASGU>

% ==== 1) Extraction des contours et infos CT ====
[domain, shapes, ext_smooth, info] = extract_slice_polys(patient_id, z_slice);

% ==== 2) Affichage CT et sauvegarde slice CT brute ====
ct_file = fullfile('Data_set', patient_id, 'ct.nii.gz');
[CT, info] = read_nifti3D(ct_file);
img = squeeze(CT(:,:,z_slice));
try
    plot_slice_ct(img, fullfile(rootOut,'slice_ct_brut.png'), info, z_slice, displayMode);
catch
    plot_slice_ct(img, fullfile(rootOut,'slice_ct_brut.png'), info, z_slice);
end

% ==== 3) Affichage slice segmentée ====
seg_dir = fullfile('Data_set', patient_id, 'segmentations');
plot_slice_segmentee(seg_dir, z_slice, ext_smooth, info, ...
    fullfile(rootOut, 'slice_segmentee.png'), displayMode);

% ==== 4) Contours extraits (avant mesh) ====
plot_contours_avant_mesh(shapes, fullfile(rootOut,'contours_avant_mesh.png'), ...
                         domain, ext_smooth, displayMode);

% ==== 5) Maillage ====
mesh_params = struct('targetSize', 5, 'minArea_mm2', 300);
meshObj = MeshBuilder.from_polys(shapes, domain, ext_smooth, z_slice, mesh_params);
meshObj.saveMesh(fullfile(rootOut, 'mesh'), z_slice);

% ==== 6) Mesh seul (VISIBLE + export) ====
fig = figure('Name','Mesh FEM seul','Color','w');
triplot(meshObj.H, meshObj.g(:,1), meshObj.g(:,2), 'k-');
axis equal tight off; title('Mesh FEM seul','FontWeight','bold');
try, hide_axes_toolbar(gca); end
apply_display_convention(gca, displayMode);
exportgraphics(fig, fullfile(rootOut, 'mesh_fem_seul.png'), 'Resolution', 300);

% ==== 7) Mesh + contours (VISIBLE + export) ====
plot_mesh_fem_contours(meshObj, fullfile(rootOut,'mesh_fem_contours.png'), displayMode, true);

% ==== 8) Électrodes (contour principal du mesh) ====
[E, el_centers, boundary_edges, boundary_nodes, contour_pts] = ... %#ok<ASGLU>
    buildElectrodesFromContour(meshObj.g, meshObj.contour, 16, 33, meshObj.H);

% ==== 9) Visualisation électrodes (VISIBLE + export) ====
fig = figure('Name','Mesh + électrodes','Color','w'); hold on; axis equal; box on;
triplot(meshObj.H, meshObj.g(:,1), meshObj.g(:,2), 'Color', [0.7 0.7 0.7]);
plot([meshObj.contour(:,1); meshObj.contour(1,1)], ...
     [meshObj.contour(:,2); meshObj.contour(1,2)], 'r-', 'LineWidth', 2);
for k = 1:numel(E)
    Ek = E{k};
    for e = 1:size(Ek,1)
        plot(meshObj.g(Ek(e,:),1), meshObj.g(Ek(e,:),2), '-', 'Color', [0.2 0.2 1], 'LineWidth', 3);
    end
end
scatter(el_centers(:,1), el_centers(:,2), 70, 'b', 'filled');
for k = 1:numel(E)
    text(el_centers(k,1), el_centers(k,2), sprintf('%d',k), ...
        'FontSize',13, 'FontWeight','bold', 'Color','b', ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
set(gca,'XTick',[],'YTick',[]);
try, hide_axes_toolbar(gca); end
apply_display_convention(gca, displayMode);
exportgraphics(fig, fullfile(rootOut, 'mesh_electrodes.png'), 'Resolution', 300);

% ==== 10) Paramètres EIT ====
eit_params = struct( ...
    'Ne',        16, ...
    'el_width',  33, ...         % mm
    'I_amp',     0.35e-3, ...    % A
    'z_contact', 0.05, ...       % Ohm·m^2
    'groups',    {meshObj.groups}, ...
    'cond',      struct('SoftTissue',0.3,'Heart',0.5,'Lung',0.15,'Trachea',0.15,'Bone',0.05,'Other',0.4), ...
    'noiseRel',  1e-3, ...
    'rngSeed',   123 ...
);

% ==== 11) Simulation EIT (forward + inverse) ====
eitSim = EITSim(meshObj.g, meshObj.H, meshObj.triGroup, meshObj.domain, eit_params, E, el_centers);
eitSim = eitSim.simulate_forward();
eitSim = eitSim.simulate_inverse([], true);  % TRUE = montrer les plots internes du solver

% ==== 12) Sigma forward (VISIBLE + export) ====
fig = figure('Name','sigma (forward)','Color','w'); hold on;
patch('Faces',eitSim.H,'Vertices',eitSim.g,'FaceVertexCData',eitSim.sigma_tri, ...
      'FaceColor','flat','EdgeColor',[0.7 0.7 0.7]);
triplot(eitSim.H, eitSim.g(:,1), eitSim.g(:,2), 'Color', [0.6 0.6 0.6]);
for k = 1:eitSim.Ne
    Ek = eitSim.E{k};
    for e = 1:size(Ek,1)
        plot(eitSim.g(Ek(e,:),1), eitSim.g(Ek(e,:),2), '-', 'Color', [0.2 0.2 1], 'LineWidth', 3);
    end
end
scatter(eitSim.el_centers(:,1), eitSim.el_centers(:,2), 70, 'b', 'filled');
axis equal tight; colormap(jet); cb=colorbar; ylabel(cb,'S/m');
title('sigma (forward)','FontWeight','bold'); set(gca,'XTick',[],'YTick',[]);
try, hide_axes_toolbar(gca); end
apply_display_convention(gca, displayMode);
exportgraphics(fig, fullfile(rootOut, 'sigma_forward_mesh_electrodes.png'), 'Resolution', 300);

% ==== 13) Sigma init (inverse) (VISIBLE + export) ====
sigma_init = ones(size(eitSim.g,1),1) * eitSim.get_cond('SoftTissue');
fig = figure('Name','sigma_{init} (inverse)','Color','w'); hold on;
patch('Faces',eitSim.H,'Vertices',eitSim.g,'FaceVertexCData',sigma_init, ...
      'FaceColor','interp','EdgeColor',[0.7 0.7 0.7]);
triplot(eitSim.H, eitSim.g(:,1), eitSim.g(:,2), 'Color', [0.7 0.7 0.7]);
for k = 1:eitSim.Ne
    Ek = eitSim.E{k};
    for e = 1:size(Ek,1)
        plot(eitSim.g(Ek(e,:),1), eitSim.g(Ek(e,:),2), '-', 'Color', [0.2 0.2 1], 'LineWidth', 3);
    end
end
scatter(eitSim.el_centers(:,1), eitSim.el_centers(:,2), 70, 'b', 'filled');
axis equal tight; colormap(jet); cb=colorbar; ylabel(cb,'S/m');
title('sigma_{init} (inverse)','FontWeight','bold'); set(gca,'XTick',[],'YTick',[]);
try, hide_axes_toolbar(gca); end
apply_display_convention(gca, displayMode);
exportgraphics(fig, fullfile(rootOut, 'sigma_init_inverse.png'), 'Resolution', 300);

% ==== 14) Plots automatiques (VISIBLE, pas d'export ici) & sauvegarde ====
eitSim.plot_all(displayMode);     % figures visibles (pas d'export)
eitSim.save_results(rootOut);

disp('Pipeline terminé.');
