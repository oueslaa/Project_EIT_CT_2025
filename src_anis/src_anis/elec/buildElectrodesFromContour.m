function [E, el_centers, boundary_edges, boundary_nodes, contour_pts] = buildElectrodesFromContour(g, contour, Ne, el_width, H)
% BUILDELECTRODESFROMCONTOUR
% -------------------------------------------------------------------------
% Place des électrodes "segmentaires" sur le contour libre (bord) d'un
% maillage 2D (PDE Toolbox).
%
% Idée
% ----
% 1) On récupère le **bord** du maillage via freeBoundary(triangulation(H,g)).
% 2) On construit la polyligne de contour ordonnée et sa longueur cumulée.
% 3) On répartit Ne centres d'électrodes à espacement régulier le long du
%    périmètre, puis on sélectionne les arêtes du bord couvertes par une
%    fenêtre d'arc de longueur `el_width` centrée sur chaque centre.
%
% Entrées
% -------
% g         : [N x 2] sommets du maillage (mm)
% H         : [M x 3] connectivité triangles (indices dans g)
% contour   : [K x 2] (mm) — non utilisé ici pour le placement (conservation
%             possible, mais le bord vient du maillage). Gardé pour compat.
% Ne        : nombre d'électrodes (entier >= 1)
% el_width  : largeur (longueur d'arc, en mm) d'une électrode
%
% Sorties
% -------
% E             : cell(Ne,1), chaque cellule contient un tableau [Le_i x 2]
%                 d'**arêtes (indices de nœuds)** appartenant à l'électrode i.
% el_centers    : [Ne x 2] centres géométriques des électrodes sur le bord (mm)
% boundary_edges: [Lb x 2] arêtes de bord (paires d'indices dans g)
% boundary_nodes: [~ x 1]  indices de nœuds présents sur le bord
% contour_pts   : [Lb x 2] points (mm) du bord dans l'ordre (première colonne
%                 de boundary_edges), pratique pour tracer la polyligne
%
% Hypothèses & remarques
% ----------------------
% - Le bord renvoyé par freeBoundary(TR) couvre le périmètre une fois.
%   Ici, on l'exploite comme une **chaîne ordonnée** d'arêtes. La ligne
%   `contour_pts = g(fb(:,1),:)` reconstruit cette polyligne.
% - `el_width` est interprété comme une **longueur d'arc**. Une électrode
%   correspond donc à toutes les arêtes du bord dont l'arc cumulé tombe
%   dans la fenêtre centrée au centre d'électrode, en gérant le **bouclage**
%   (modulo périmètre).
% - Si `el_width` est trop grand, des électrodes peuvent se chevaucher
%   (c'est permis : chaque cellule E{k} est indépendante).
%
% Exemple
% -------
% TR = triangulation(H,g); figure; triplot(TR); axis equal
% [E, C, be, bn, cps] = buildElectrodesFromContour(g, [], 16, 12, H);
% hold on; plot(cps(:,1),cps(:,2),'k.-'); plot(C(:,1),C(:,2),'ro');
% % Visualiser la k-ième électrode:
% k=1; ed = E{k}; for t=1:size(ed,1), plot(g(ed(t,:),1), g(ed(t,:),2), 'LineWidth',3); end
%
% -------------------------------------------------------------------------

% --- 1) Bord du maillage --------------------------------------------------
TR = triangulation(H, g);
fb = freeBoundary(TR);            % [Lb x 2] : arêtes de bord (chaîne ordonnée)
boundary_edges = fb;              % pour retour direct
boundary_nodes = unique(fb(:));   % tous les nœuds impliqués dans le bord
contour_pts    = g(fb(:,1), :);   % polyligne ordonnée (points correspondant à fb(:,1))

% --- 2) Longueurs d'arc cumulées le long du bord --------------------------
arc_len = [0; cumsum(vecnorm(diff(contour_pts),2,2))];
perim   = arc_len(end);           % périmètre total (mm)
spacing = perim / Ne;             % espacement régulier des centres d'électrodes

% Pré-allocation des sorties principales
E = cell(Ne,1);
el_centers = zeros(Ne,2);

% --- 3) Boucle sur les électrodes -----------------------------------------
for k = 1:Ne
    % 3.a) Position d'arc du centre (mm depuis l'origine de la polyligne)
    center_pos = (k-1)*spacing + spacing/2;

    % Trouver le segment [i -> j] de la polyligne qui contient center_pos
    i = find(arc_len <= center_pos, 1, 'last');
    j = i + 1; if j>numel(arc_len), j=2; end

    % Interpolation linéaire pour le centre en coordonnées (x,y)
    t = (center_pos - arc_len(i)) / max(arc_len(j)-arc_len(i), eps);
    el_centers(k,:) = contour_pts(i,:) + t*(contour_pts(j,:)-contour_pts(i,:));

    % 3.b) Fenêtre d'arc pour l'électrode k : [center - w/2, center + w/2]
    % Gestion du **wrap-around** avec modulo périmètre.
    a0 = mod(center_pos - el_width/2, perim);
    a1 = mod(center_pos + el_width/2, perim);

    % in(s) indique si une position d'arc s (en mm) tombe dans la fenêtre.
    % Cas 1 (a0<a1) : intervalle standard    -> a0 <= s <= a1
    % Cas 2 (a0>=a1) : intervalle "qui déborde" -> s >= a0 OU s <= a1
    in = @(s) (a0<a1).*(s>=a0 & s<=a1) + (a0>=a1).*(s>=a0 | s<=a1);

    % On sélectionne les **arêtes** du bord dont AU MOINS une extrémité
    % de la polyligne tombe dans la fenêtre (approx. robuste).
    mask = in(arc_len(1:end-1)) | in(arc_len(2:end));

    % fb(mask,:) renvoie les paires d'indices de nœuds pour l'électrode k
    sel  = fb(mask,:);
    E{k} = sel;  % arêtes (paires d'indices) du bord pour l'électrode k
end
end
