function plot_mesh_fem_contours(meshObj, outpath, mode, show)
% PLOT_MESH_FEM_CONTOURS  Mesh coloré par groupes + contours organes
% Usage:
%   plot_mesh_fem_contours(meshObj, outpath, mode, show)
%   - mode (opt) : 'radiological' (defaut) ou 'neurological'
%   - show (opt) : true (defaut) -> figure visible; false -> invisible (batch)

if nargin < 3 || isempty(mode), mode = 'radiological'; end
if nargin < 4 || isempty(show), show = true; end

vis = 'on'; if ~show, vis = 'off'; end
fig = figure('Visible', vis, 'Color','w', ...
             'Name','Mesh FEM + contours','NumberTitle','off');
hold on; axis equal off;

% Palette (6 groupes)
cmap = [
    0.4 0.7 1.0;   % Soft Tissue
    1.0 0.6 0.3;   % Heart
    1.0 0.9 0.4;   % Lung
    0.4 1.0 0.7;   % Trachea
    0.5 0.85 1.0;  % Bone
    0.5 0.7 0.5    % Other
];

% Triangles par groupe
for k = 1:numel(meshObj.groups)
    idx = (meshObj.triGroup == k);
    if ~any(idx), continue; end
    tris = meshObj.H(idx,:);
    patch('Faces', tris, 'Vertices', meshObj.g, ...
          'FaceColor', cmap(k,:), 'FaceAlpha', 0.85, ...
          'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.6);
end

% Contour extérieur
P = meshObj.contour;
if ~isequal(P(1,:), P(end,:)), P = [P; P(1,:)]; end
plot(P(:,1), P(:,2), 'k-', 'LineWidth', 3);

% Contours d'organes (sauf soft tissue)
for i = 2:numel(meshObj.groups)
    cc  = cmap(i,:);
    shp = meshObj.shapes.(meshObj.groups{i});
    if shp.NumRegions == 0, continue; end
    regs = regions(shp);
    for r = 1:numel(regs)
        V = regs(r).Vertices;
        if ~isequal(V(1,:), V(end,:)), V = [V; V(1,:)]; end
        plot(V(:,1), V(:,2), '-', 'Color', cc, 'LineWidth', 2);
    end
end

title('Mesh FEM + contours (groupes)', 'Interpreter','none','FontWeight','bold','FontSize',16);
set(gca,'XTick',[],'YTick',[]); axis tight;

% Orientation + toolbar
try, hide_axes_toolbar(gca); end %#ok<TRYNC>
apply_display_convention(gca, mode);

% Export si demandé
if nargin >= 2 && ~isempty(outpath)
    exportgraphics(fig, outpath, 'Resolution', 300);
end

% Si visible, on force l'affichage
if strcmpi(vis,'on'), drawnow; else, close(fig); end
end
