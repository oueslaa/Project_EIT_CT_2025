function [E, el_centers, boundary_edges, boundary_nodes, contour_pts] = buildElectrodesFromContour(g, contour, Ne, el_width, H)
% Place les électrodes sur le **contour libre ordonné** du maillage
% g [N x 2] en mm, H [M x 3], Ne (#électrodes), el_width en mm

TR = triangulation(H, g);
fb = freeBoundary(TR);            % Nx2 (chaîne ordonnée de noeuds)
boundary_edges = fb;              % déjà des arêtes de bord
boundary_nodes = unique(fb(:));
contour_pts    = g(fb(:,1), :);  % ordonné

arc_len = [0; cumsum(vecnorm(diff(contour_pts),2,2))];
perim   = arc_len(end);
spacing = perim / Ne;

E = cell(Ne,1);
el_centers = zeros(Ne,2);

for k = 1:Ne
    center_pos = (k-1)*spacing + spacing/2;
    i = find(arc_len <= center_pos, 1, 'last');
    j = i + 1; if j>numel(arc_len), j=2; end
    t = (center_pos - arc_len(i)) / max(arc_len(j)-arc_len(i), eps);
    el_centers(k,:) = contour_pts(i,:) + t*(contour_pts(j,:)-contour_pts(i,:));

    a0 = mod(center_pos - el_width/2, perim);
    a1 = mod(center_pos + el_width/2, perim);
    in = @(s) (a0<a1).*(s>=a0 & s<=a1) + (a0>=a1).*(s>=a0 | s<=a1);
    mask = in(arc_len(1:end-1)) | in(arc_len(2:end));
    sel  = fb(mask,:);
    E{k} = sel;  % arêtes (paires d'indices) sur le bord
end
end