function plot_contours_avant_mesh(shapes, outpng, domain, ext_smooth, mode)
% PLOT_CONTOURS_AVANT_MESH  Affiche les contours par organes + contour externe
% Usage:
%   plot_contours_avant_mesh(shapes, outpng, domain, ext_smooth, mode)
%   - mode (optionnel): 'radiological' (defaut) ou 'neurological'
%
% Dependances: apply_display_convention(ax,mode) + (optionnel) hide_axes_toolbar(ax)

if nargin < 5 || isempty(mode), mode = 'radiological'; end

fig = figure('Color','w'); hold on; axis equal; axis off;

% Couleurs par organe
colors = struct('Heart',[0.15 0.5 1], 'Lung',[1 0.4 0], 'Trachea',[0.5 0 1], ...
                'Bone',[0.3 0.9 0.3], 'Other',[0.5 0.7 0.5]);

fields = fieldnames(shapes);
legend_handles = [];
legend_labels  = {};

for i = 1:numel(fields)
    organ = fields{i};
    shp   = shapes.(organ);
    if shp.NumRegions == 0, continue; end
    cc = colors.Other; if isfield(colors, organ), cc = colors.(organ); end
    regs = regions(shp);
    for r = 1:numel(regs)
        P = regs(r).Vertices;
        if ~isequal(P(1,:), P(end,:)), P = [P; P(1,:)]; end
        h = plot(P(:,1), P(:,2), '-', 'LineWidth', 2, 'Color', cc);
        if r == 1
            legend_handles(end+1) = h; %#ok<AGROW>
            legend_labels{end+1}  = organ; %#ok<AGROW>
        end
    end
end

% Contour externe
if exist('ext_smooth','var') && ~isempty(ext_smooth)
    P = ext_smooth;
    if ~isequal(P(1,:), P(end,:)), P = [P; P(1,:)]; end
    h_ext = plot(P(:,1), P(:,2), 'k-', 'LineWidth', 4);
    legend_handles(end+1) = h_ext; legend_labels{end+1} = 'Exterieur';
elseif exist('domain','var') && ~isempty(domain) && isa(domain,'polyshape')
    [x,y] = boundary(domain); P = [x(:), y(:)];
    if ~isequal(P(1,:), P(end,:)), P = [P; P(1,:)]; end
    h_ext = plot(P(:,1), P(:,2), 'k-', 'LineWidth', 4);
    legend_handles(end+1) = h_ext; legend_labels{end+1} = 'Exterieur';
end

title('Contours extraits (avant mesh)', 'FontWeight','bold','FontSize',16);
legend(legend_handles, legend_labels, 'Location','bestoutside');

% Orientation + toolbar
try, hide_axes_toolbar(gca); end %#ok<TRYNC>
apply_display_convention(gca, mode);

exportgraphics(fig, outpng, 'Resolution', 220);
end
