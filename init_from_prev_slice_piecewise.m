function [sigma0_nod, info] = init_from_prev_slice_piecewise(meshPrevFile, packPrevFile, meshThisFile, opts)
% σ0 piecewise-constant: transfert des "couleurs" (valeurs discrètes)
% de la slice PRECEDENTE vers la slice COURANTE (contour/électrodes de la cible).
%
% Entrées
%   meshPrevFile : Outputs/.../slice_302/mesh/mesh_slice302.mat  (contour 302)
%   packPrevFile : Outputs/.../slice_302/eit_pack.mat            (g,H,E,sigma_tri 302)
%   meshThisFile : Outputs/.../slice_301/mesh/mesh_slice301.mat  (contour 301)
%   opts.doWarp  : (bool) aligne les contours 301→302 (similarité) [def=true]
%   opts.clip    : [min max] bornes sigma                         [def=[0 inf]]
%
% Sorties
%   sigma0_nod : (Nnod x 1) nodal, piecewise-constant par organe
%   info       : struct (niveaux, labels, correspondances, etc.)

if nargin<4 || ~isstruct(opts), opts = struct; end
doWarp = true;  if isfield(opts,'doWarp'), doWarp = opts.doWarp; end
clip   = [0, inf]; if isfield(opts,'clip'), clip = opts.clip; end
quantTol = 1e-3;  % tolérance pour regrouper les niveaux constants

% --------- charge fichiers ---------
MP = load(meshPrevFile);       % pour les contours: MP.contour
PP = load(packPrevFile);       % gP,H P, sigma_tri (302)
MT = load(meshThisFile);       % pour le contour cible: MT.contour

% maillages
gP = PP.g;   HP = PP.H;   % prev (302)
gT = MT.g;   HT = MT.H;   % this (301)

% --------- niveaux (couleurs) de la slice 302 ---------
assert(isfield(PP,'sigma_tri'), 'packPrevFile doit contenir sigma_tri (valeur tri).');
sigp_tri = PP.sigma_tri(:);
[levels, labP] = levels_from_sigma(sigp_tri, quantTol);  % niveaux distincts + labels tri (302)

% --------- aligne 301 → 302 (optionnel) ---------
gT_in_prev = gT;
if doWarp && isfield(MP,'contour') && isfield(MT,'contour') ...
          && ~isempty(MP.contour) && ~isempty(MT.contour)
    Cprev = MP.contour;
    Cthis = MT.contour;
    gT_in_prev = map_this_to_prev(gT, Cprev, Cthis);      % 301 -> repère 302
end

% --------- affecte un label à CHAQUE triangle de 301 ---------
cT = tri_centroids(gT, HT);               % centroïdes tri 301 (dans repère 302)
TR = triangulation(HP, gP);               % triangulation source (302)
tid = pointLocation(TR, cT);              % tri 302 contenant chaque centroïde 301 (NaN si hors)

labT = zeros(size(HT,1),1);               % labels pour tri 301
in    = ~isnan(tid);
labT(in) = labP(tid(in));                 % direct quand le centroïde tombe dans un tri 302

% pour les centroïdes hors-domaine: plus proche centroïde de 302
if any(~in)
    cP = tri_centroids(gP, HP);
    miss = ~in;
    idxNN = knn_fallback(cP, cT(miss,:)); % index tri 302 le plus proche
    labT(miss) = labP(idxNN);
end

% --------- convertit labels TRI -> valeurs NODES par vote majoritaire ---------
sigma0_nod = labels_to_nodes(labT, HT, levels);

% --------- bornes ---------
sigma0_nod = max(sigma0_nod, clip(1));
sigma0_nod = min(sigma0_nod, clip(2));

% --------- info ---------
info = struct('levels',levels,'labels_prev_tri',labP,'labels_this_tri',labT, ...
              'centroids_this',cT,'tid_in_prev',tid);
end

% ==================== helpers locaux ====================

function [levels, labels] = levels_from_sigma(sigTri, tol)
% regroupe des valeurs constantes (par organes) avec une tolérance
if nargin<2 || isempty(tol), tol = 1e-3; end
r = round(sigTri(:)/tol)*tol;           % quantification robuste
[levels, ~, labels] = unique(r);       % labels ∈ {1..K}
levels = sort(levels(:));
% remet labels en cohérence avec l'ordre trié des levels
[~,loc] = ismember(r, levels);
labels = loc(:);
end

function Xp = map_this_to_prev(X, Cprev, Cthis)
% mappe des points de la géométrie THIS (301) vers le repère PREV (302)
try
    % fitgeotrans(moving, fixed): moving=Cthis -> fixed=Cprev
    tform = fitgeotrans(Cthis, Cprev, 'similarity');
    Xp = transformPointsForward(tform, X);
catch
    % fallback: Procrustes (transforme Cthis -> Cprev)
    [~,~,T] = procrustes(Cprev, Cthis, 'scaling', true, 'reflection', false);
    % Z = T.b * Y * T.T + T.c ; on applique aux X (Y-espace)
    Xp = T.b * X * T.T + T.c(1,:);
end
end

function C = tri_centroids(V, T)
C = (V(T(:,1),:) + V(T(:,2),:) + V(T(:,3),:)) / 3;
end

function idx = knn_fallback(P, Q)
% renvoie pour chaque point de Q l'indice du point le plus proche dans P
try
    idx = knnsearch(P, Q);
catch
    % simple (O(NM)) mais robuste
    nQ = size(Q,1);
    idx = zeros(nQ,1);
    for i = 1:nQ
        d = sum((P - Q(i,:)).^2, 2);
        [~, idx(i)] = min(d);
    end
end
end

function sigma_nod = labels_to_nodes(labTri, T, levels)
% assigne à chaque NOEUD la valeur du niveau majoritaire parmi ses triangles incidents
nNode = max(T(:));
nTri  = size(T,1);
I = [T(:,1); T(:,2); T(:,3)];
J = [(1:nTri)'; (1:nTri)'; (1:nTri)'];
cellTri = accumarray(I, J, [nNode 1], @(x){x});
sigma_nod = zeros(nNode,1);
for i=1:nNode
    tr = cellTri{i};
    if isempty(tr), continue; end
    labs = labTri(tr);
    % mode (vote majoritaire) — si ex-aequo, 'mode' prend le plus petit label
    lab = mode(labs);
    sigma_nod(i) = levels(lab);
end
end
