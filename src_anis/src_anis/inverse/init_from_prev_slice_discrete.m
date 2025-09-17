function [sigma0_inv, info] = init_from_prev_slice_discrete(meshPrevFile, packPrevFile, invMeshFile, packThisFile, opts)
% INIT_FROM_PREV_SLICE_DISCRETE
% -------------------------------------------------------------------------
% Initialise la conductivité NODALE de la slice cible (maillage INVERSE)
% en transplantant les organes (polyshapes) de la slice précédente et en
% les plaçant dans le contour extérieur cible **comme le MeshBuilder** :
%   1) Warp (TPS ou similarité) des organes vers le contour cible,
%   2) Lissage & simplification des polygones, clip dans le domaine cible,
%   3) Exclusivité stricte des organes (ordre de priorité),
%   4) Étiquetage nodal + snapping bande + nettoyage d'îlots,
%   5) Mapping labels -> σ via params.cond de la cible.
%
% Entrées
%   meshPrevFile : .mat slice précédente (g,H,contour, shapes{group})
%   packPrevFile : réservé (non requis ici)
%   invMeshFile  : .mat maillage INVERSE cible (g,H,contour)
%   packThisFile : .mat pack FORWARD cible (params.cond)
%   opts         : struct optionnelle (voir `DEFAULTS` ci-dessous)
%
% Sorties
%   sigma0_inv : [NnI x 1] σ nodale sur le maillage INVERSE
%   info       : diagnostics (labels, levels, warp_used, shapes_warped, ...)
%
% ©2025

% ---------- Chargements ----------
Mprev = load(meshPrevFile);   %#ok<NASGU>
Minv  = load(invMeshFile);
Stgt  = load(packThisFile);

assert(isfield(Minv,'g') && isfield(Minv,'H') && isfield(Minv,'contour'), ...
    'invMeshFile doit contenir g,H,contour.');
assert(isfield(Stgt,'params') && isfield(Stgt.params,'cond'), ...
    'packThisFile doit contenir params.cond.');

gI_mm = Minv.g; HI = Minv.H; NnI = size(gI_mm,1);

% ---------- DEFAULTS ----------
if nargin < 5, opts = struct(); end
opts.clip                 = local_def(opts,'clip', [1e-3, 2.0]);
opts.doWarp               = local_def(opts,'doWarp', true);

opts.init                 = local_def(opts,'init', struct());
opts.init.keep_groups     = local_def(opts.init,'keep_groups', ...
    {'Heart','Lung','Trachea','Bone','Cartilage','Blood','LiverSpleenPancreas','Kidney','Other'});
opts.init.knn_k           = local_def(opts.init,'knn_k', 11);
opts.init.clean_iters     = local_def(opts.init,'clean_iters', 2);
opts.init.border_smooth_mm= local_def(opts.init,'border_smooth_mm', 1.2);
opts.init.min_area_mm2    = local_def(opts.init,'min_area_mm2', 150);
opts.init.contour_simplify_tol = local_def(opts.init,'contour_simplify_tol', 0.6);
opts.init.tau_uncertainty = local_def(opts.init,'tau_uncertainty', 0.90);

opts.tps                  = local_def(opts,'tps', struct());
opts.tps.enable           = local_def(opts.tps,'enable', true);
opts.tps.lambda           = local_def(opts.tps,'lambda', 1e-3);

opts.snap                 = local_def(opts,'snap', struct());
opts.snap.band_mm         = local_def(opts.snap,'band_mm', 2.5);
opts.snap.island_min_area_mm2 = local_def(opts.snap,'island_min_area_mm2', 150);

% ---------- Groupes / domaine cible ----------
full_groups  = {'Heart','Lung','Trachea','Bone','Cartilage','Blood','LiverSpleenPancreas','Kidney','Other'};
keep_groups  = opts.init.keep_groups;
groups_order = full_groups(ismember(full_groups, keep_groups));

Pext = Minv.contour; 
if ~isequal(Pext(1,:),Pext(end,:)), Pext = [Pext; Pext(1,:)]; end
domain_target = polyshape(Pext(:,1), Pext(:,2), 'Simplify', true);

% ---------- 1) Shapes d'organes (slice précédente) ----------
shPrev = local_get_shapes_from_prev(load(meshPrevFile), groups_order);

% ---------- 2) Warp vers le contour cible ----------
if opts.doWarp
    [shWarp, warp_used] = local_warp_shapes(shPrev, domain_target, opts, load(meshPrevFile), Minv);
else
    shWarp = shPrev; warp_used = false;
end

% ---------- 2a) Placement "façon MeshBuilder" ----------
shUse = local_place_organs_meshbuilder( ...
    shWarp, Minv.contour, ...
    struct('smooth_mm',  opts.init.border_smooth_mm, ...
           'minArea',    opts.init.min_area_mm2, ...
           'simpl_tol',  opts.init.contour_simplify_tol), ...
    groups_order);


% ---------- 3) Étiquetage nodal + snapping + nettoyage ----------
XY = gI_mm; % mm
labels = local_label_nodes_by_shapes(XY, shUse, groups_order, opts.init.knn_k, opts.init.clean_iters);
labels = local_snap_labels_by_distance(labels, shUse, XY, opts.snap.band_mm, groups_order);
labels = local_remove_islands(labels, gI_mm, HI, opts.snap.island_min_area_mm2);

% ---------- 4) labels -> σ nodale (params.cond cible) ----------
sigma0_inv = local_labels_to_sigma(labels, groups_order, Stgt.params, opts.clip);
sigma0_inv = sigma0_inv(:);

% ---------- Infos sortie ----------
info = struct();
info.lab_nod     = labels(:);
info.lab_tri     = []; % si besoin plus tard
info.levels      = local_levels_from_cond(Stgt.params.cond, groups_order);
info.warp_used   = warp_used;
info.init_mode   = 'morpho';
info.shapes_used = shUse;

if numel(sigma0_inv) ~= NnI
    error('init_from_prev_slice_discrete: taille σ0 (%d) != NnI (%d).', numel(sigma0_inv), NnI);
end
end

%% ===================== UTILITAIRES — SHAPES & WARP ======================

function sh = local_get_shapes_from_prev(Mprev, groups_order)
% Retourne un struct de polyshape pour chaque groupe demandé (sans trous).
    sh = struct();
    hasShapes = isfield(Mprev,'shapes');
    for k = 1:numel(groups_order)
        gk = groups_order{k};
        if hasShapes && isfield(Mprev.shapes, gk)
            P = Mprev.shapes.(gk);
            P = local_rmholes(P);
            sh.(gk) = simplify(P,'KeepCollinear',true);
        else
            warning('init: shapes.%s absent — utilisé vide.', gk);
            sh.(gk) = polyshape();
        end
    end
end

function [shW, warp_used] = local_warp_shapes(shPrev, domain_target, opts, Mprev, Minv)
% Aligne la slice précédente sur la cible : TPS (si enable) sinon similarité.
    warp_used = false;

    Cprev = local_chain_from_tri(Mprev.g, Mprev.H);
    if isempty(Cprev) && isfield(Mprev,'contour'), Cprev = Mprev.contour; end
    if isempty(Cprev)
        warning('Contour précédent indisponible — skip warp.');
        shW = shPrev; return;
    end
    if ~isequal(Cprev(1,:), Cprev(end,:)), Cprev = [Cprev; Cprev(1,:)]; end
    Ctgt  = Minv.contour; if ~isequal(Ctgt(1,:), Ctgt(end,:)), Ctgt=[Ctgt;Ctgt(1,:)]; end

    Mmatch = min(max(150, size(Cprev,1)), 800);
    Psrc = local_resample_polyline(Cprev, Mmatch);
    Ptgt = local_resample_polyline(Ctgt,  Mmatch);

    try
        if opts.tps.enable
            modelTPS = local_fit_tps(Psrc, Ptgt, opts.tps.lambda);
            shW      = local_warp_tps_shapes(shPrev, modelTPS);
        else
            Taff     = local_fit_similarity(Psrc, Ptgt);
            shW      = local_affine_shapes(shPrev, Taff);
        end
        warp_used = true;
        % clip préventif (sécurité)
        flds = fieldnames(shW);
        for i=1:numel(flds)
            if shW.(flds{i}).NumRegions==0, continue; end
            shW.(flds{i}) = simplify(intersect(shW.(flds{i}), domain_target), 'KeepCollinear',true);
        end
    catch ME
        warning('Warp a échoué (%s) — shapes non-warpés.', ME.message);
        shW = shPrev;
    end
end

function P2 = local_rmholes(P1)
% Enlève les trous **sans** 'boundaries' (compatible anciennes versions).
    if P1.NumRegions == 0
        P2 = P1; return;
    end
    regs = regions(P1);
    P2 = polyshape();
    for r = 1:numel(regs)
        [x,y] = local_region_outer_xy(regs(r));
        if numel(x) >= 3
            Q = polyshape(x,y,'Simplify',true);
            P2 = union(P2, Q);
        end
    end
    P2 = simplify(P2,'KeepCollinear',true);
end

function [x,y] = local_region_outer_xy(R)
% Renvoie la boucle **externe** d'un polyshape région R via Vertices.
    V = R.Vertices;
    if isempty(V)
        x = []; y = []; return;
    end
    isn = any(isnan(V),2);
    if any(isn)
        k = find(isn,1,'first')-1;
        if isempty(k) || k < 3
            x = []; y = []; return;
        end
        V = V(1:k,:);
    end
    if size(V,1) >= 3
        if ~isequal(V(1,:), V(end,:)), V = [V; V(1,:)]; end
        x = V(:,1); y = V(:,2);
    else
        x = []; y = [];
    end
end

function shI = local_place_organs_meshbuilder(shIn, ext_target, p, groups_order)
% Place les organes "comme le mesher" :
%  - lissage (polybuffer closing/opening) + simplification,
%  - clip dans le domaine,
%  - exclusivité stricte des organes (ordre groups_order).
    if nargin<3 || isempty(p), p = struct(); end
    if ~isfield(p,'smooth_mm'),  p.smooth_mm  = 1.2; end
    if ~isfield(p,'minArea'),    p.minArea    = 150; end
    if ~isfield(p,'simpl_tol'),  p.simpl_tol  = 0.5; end

    % Domaine (contour ouvert -> polyshape)
    P = ext_target;
    if ~isequal(P(1,:),P(end,:)), P=[P;P(1,:)]; end
    domain_target = polyshape(P(:,1), P(:,2), 'Simplify', true);

    % 1) Lissage + simplification
    shClean = struct();
    flds = fieldnames(shIn);
    for i=1:numel(flds)
        Pi = shIn.(flds{i});
        if Pi.NumRegions==0, shClean.(flds{i}) = polyshape(); continue; end
        try % closing + opening
            Pi = polybuffer(polybuffer(Pi, +p.smooth_mm, 'JointType','round'), -p.smooth_mm, 'JointType','round');
            Pi = polybuffer(polybuffer(Pi, -p.smooth_mm, 'JointType','round'), +p.smooth_mm, 'JointType','round');
        catch
            % polybuffer pas dispo => on garde Pi tel quel
        end
        regs = regions(Pi); Qall = polyshape();
        for r=1:numel(regs)
            [x,y] = local_region_outer_xy(regs(r));
            if numel(x)<3, continue; end
            xy = local_simplify_ring([x y], p.simpl_tol);
            xy = local_clean_ring(xy, 1e-6);
            Q  = polyshape(xy(:,1), xy(:,2), 'Simplify', true);
            if area(Q) >= p.minArea
                Qall = union(Qall, Q);
            end
        end
        shClean.(flds{i}) = simplify(Qall,'KeepCollinear',true);
    end

    % 2) Clip domaine
    shClip = struct();
    for i=1:numel(flds)
        Q = shClean.(flds{i});
        if Q.NumRegions==0, shClip.(flds{i}) = polyshape(); continue; end
        shClip.(flds{i}) = simplify(intersect(Q, domain_target), 'KeepCollinear', true);
    end

    % 3) Exclusivité par priorité
    shI = local_make_shapes_exclusive(shClip, groups_order);
end

function Sx = local_make_shapes_exclusive(S, order)
% Règle d'exclusivité : on remplit dans l'ordre donné et on "creuse" les suivants.
    Sx = struct();
    Ufilled = polyshape();   % déjà occupé
    for j = 1:numel(order)
        g = order{j};
        if ~isfield(S, g) || S.(g).NumRegions==0
            Sx.(g) = polyshape(); continue;
        end
        P = subtract(S.(g), Ufilled);
        P = simplify(P, 'KeepCollinear', true);
        Sx.(g) = P;
        Ufilled = union(Ufilled, P);
    end
end

function T = local_fit_similarity(P, Q)
% Ajuste x -> s*R*x + t (Procrustes sans réflexion).
    assert(size(P,2)==2 && size(Q,2)==2 && size(P,1)==size(Q,1), 'fit_similarity: tailles incohérentes');
    Pc = mean(P,1); Qc = mean(Q,1);
    X = P - Pc; Y = Q - Qc;
    [U,~,V] = svd(X.'*Y);
    R = U*V.'; if det(R)<0, U(:,end)=-U(:,end); R=U*V.'; end
    s = trace((R.'*(X.'*Y))) / max(trace(X.'*X), eps);
    if ~isfinite(s) || s<=0, s=1; end
    t = Qc.' - s*R*Pc.';
    T = eye(3); T(1:2,1:2) = s*R; T(1:2,3) = t;
end

function shW = local_affine_shapes(sh, T)
% Applique une affine homogène 3x3 aux shapes (boucles externes).
    flds = fieldnames(sh);
    shW = struct();
    for i=1:numel(flds)
        P = sh.(flds{i});
        if P.NumRegions==0, shW.(flds{i}) = P; continue; end
        regs = regions(P);
        Qall = polyshape();
        for r=1:numel(regs)
            [x,y] = local_region_outer_xy(regs(r));
            if numel(x) < 3, continue; end
            Vh = [x(:), y(:), ones(numel(x),1)] * T.';
            Vh = local_clean_ring(Vh(:,1:2), 1e-6);      % <<< nettoie avant polyshape
            Q  = polyshape(Vh(:,1), Vh(:,2), 'Simplify', true);
            Qall = union(Qall, Q);
        end
        shW.(flds{i}) = local_rmholes(simplify(Qall, 'KeepCollinear', true));
    end
end

function C = local_chain_from_tri(g, H)
% Chaîne de bord ordonnée à partir du freeBoundary (approx).
    try
        TR = triangulation(H, g);
        fb = freeBoundary(TR);
        if isempty(fb), C = []; return; end
        C = g(fb(:,1),:);
    catch
        C = [];
    end
end

function R = local_resample_polyline(P, M)
% Échantillonne M points sur une polyligne FERMÉE P (Nx2).
    if ~isequal(P(1,:), P(end,:)), P=[P; P(1,:)]; end
    s = [0; cumsum(vecnorm(diff(P),2,2))];
    L = s(end); if L<=0, R = repmat(P(1,:),M,1); return; end
    su = linspace(0, L, M).';
    R = interp1(s, P, su, 'linear');
end

% ------------------ TPS (thin plate splines) -----------------------------

function model = local_fit_tps(P, Q, lambda)
% P,Q : Nx2. Solve (K+λI)w + P1*a = Q; P1=[P 1]; contraintes P1' w = 0.
    N = size(P,1);
    K = local_tps_K(P,P);
    K(1:N+1:end) = K(1:N+1:end) + lambda;
    P1 = [P, ones(N,1)];
    A  = [K, P1; P1.', zeros(3,3)];
    Bx = [Q(:,1); 0;0;0];
    By = [Q(:,2); 0;0;0];
    solx = A \ Bx;  soly = A \ By;
    model = struct('P',P, 'wX',solx(1:N), 'ax',solx(N+1:N+3), ...
                        'wY',soly(1:N), 'ay',soly(N+1:N+3));
end

function K = local_tps_K(X, Y)
% K_ij = r^2 log r  (r=||Xi - Yj||2), K_ii = 0.
    NX = size(X,1); NY = size(Y,1);
    K = zeros(NX,NY);
    for i=1:NX
        d = hypot(X(i,1)-Y(:,1), X(i,2)-Y(:,2));
        d(d==0) = 1; % évite log(0)
        v = (d.^2).*log(d);
        v(~isfinite(v)) = 0;
        K(i,:) = v;
    end
    if NX==NY
        ii = 1:NX; K(sub2ind([NX,NY],ii,ii)) = 0;
    end
end

function Y = local_tps_apply(model, X)
% Warp des points X (Mx2) -> Y (Mx2).
    M = size(X,1);
    K = local_tps_K(X, model.P);
    X1 = [X, ones(M,1)];
    Yx = K*model.wX + X1*model.ax;
    Yy = K*model.wY + X1*model.ay;
    Y  = [Yx, Yy];
end

function shW = local_warp_tps_shapes(sh, model)
% TPS appliqué aux boucles externes de chaque polyshape.
    flds = fieldnames(sh); shW = struct();
    for i=1:numel(flds)
        P = sh.(flds{i});
        if P.NumRegions==0, shW.(flds{i}) = P; continue; end
        regs = regions(P);
        Qall = polyshape();
        for r=1:numel(regs)
            [x,y] = local_region_outer_xy(regs(r));
            if numel(x) < 3, continue; end
            xy2 = local_tps_apply(model, [x, y]);
            xy2 = local_clean_ring(xy2, 1e-6);             % <<< crucial pour éviter warnings
            Q = polyshape(xy2(:,1), xy2(:,2), 'Simplify', true);
            Qall = union(Qall, Q);
        end
        shW.(flds{i}) = local_rmholes( simplify(Qall,'KeepCollinear',true) );
    end
end

%% ===================== UTILITAIRES — LABELLISATION ======================

function labels = local_label_nodes_by_shapes(XY, shUse, groups_order, knn_k, clean_iters)
% Étiquette nodale par priorité des groupes. Nettoyage local minimal.
    N = size(XY,1);
    labels = zeros(N,1); % 0=SoftTissue
    for j = 1:numel(groups_order)
        gj = groups_order{j};
        shp = shUse.(gj);
        if shp.NumRegions==0, continue; end
        IN = isinterior(shp, XY(:,1), XY(:,2));
        labels(IN) = j;
    end
    % (placeholder) Nettoyage supplémentaire possible ici si besoin
    if clean_iters > 0 && knn_k > 0
        % à brancher si nécessaire (majorité KNN / voisinage nodal)
    end
end

function labels = local_snap_labels_by_distance(labels, shUse, XY, band_mm, groups_order)
% Dans une bande autour des frontières (band_mm), on assigne le groupe le plus proche.
    if band_mm <= 0, return; end
    m = size(XY,1);
    G = inf(m, numel(groups_order));
    for j = 1:numel(groups_order)
        gj = groups_order{j}; shp = shUse.(gj);
        if shp.NumRegions==0, continue; end
        [dx,dy] = local_distance_to_polyshape(shp, XY);
        G(:,j) = hypot(dx,dy);
    end
    [dmin, jmin] = min(G, [], 2, 'omitnan');
    inBand = dmin <= band_mm;
    labels(inBand) = jmin(inBand);
end

function [dx,dy] = local_distance_to_polyshape(P, XY)
% Distance au bord (non signée) par projection segmentaire des boucles externes.
    regs = regions(P);
    B = [];
    for r=1:numel(regs)
        [x,y] = local_region_outer_xy(regs(r));
        if numel(x) >= 2
            B = [B; [x(:), y(:)]]; %#ok<AGROW>
        end
    end
    if isempty(B)
        dx = zeros(size(XY,1),1); dy = dx; return;
    end
    E = [ (1:size(B,1)-1).', (2:size(B,1)).' ];
    dx = zeros(size(XY,1),1); dy = dx; d2 = inf(size(XY,1),1);
    for e = 1:size(E,1)
        a = B(E(e,1),:); b = B(E(e,2),:);
        ab = b-a; L2 = sum(ab.^2); if L2==0, continue; end
        t = ((XY-a).*ab) * [1;1] / L2; t = max(0, min(1, t));
        proj = a + t.*ab;
        v = XY - proj;
        d2e = sum(v.^2,2);
        pick = d2e < d2;
        dx(pick) = v(pick,1); dy(pick) = v(pick,2); d2(pick) = d2e(pick);
    end
end

function labels = local_remove_islands(labels, g, H, min_area_mm2)
% Supprime les petits îlots par label (aire triangulaire).
    if min_area_mm2 <= 0, return; end
    Mt = size(H,1);
    % aire des triangles
    Atri = zeros(Mt,1);
    for t=1:Mt
        v = g(H(t,:),:);
        Atri(t) = 0.5*abs(det([v(2,:)-v(1,:); v(3,:)-v(1,:)]));
    end
    % label des triangles = majorité des noeuds
    labTri = zeros(Mt,1);
    for t=1:Mt
        L = labels(H(t,:));
        u = unique(L); cnt = arrayfun(@(x) sum(L==x), u);
        [~,ix] = max(cnt); labTri(t) = u(ix);
    end
    % adjacency triangle-triangle (partage d'une arête)
    AdjT = local_triangle_adjacency(H);
    labs = setdiff(unique(labTri), 0).';
    for lab = labs
        S = find(labTri==lab);
        if isempty(S), continue; end
        Sub = AdjT(S,S);
        Gsub = graph(Sub,'OmitSelfLoops');
        comp = conncomp(Gsub);
        for c = unique(comp)
            tri_c = S(comp==c);
            area_c = sum(Atri(tri_c));
            if area_c < min_area_mm2
                nodes = unique(H(tri_c,:));
                labels(nodes) = 0; % remet en SoftTissue
            end
        end
    end
end

function AdjT = local_triangle_adjacency(H)
% Matrice d'adjacence des triangles qui partagent une arête.
    Mt = size(H,1);
    edges = [H(:,[1 2]); H(:,[2 3]); H(:,[3 1])];
    triid = [(1:Mt)'; (1:Mt)'; (1:Mt)'];
    edges = sort(edges,2);
    [~,~,ic] = unique(edges, 'rows');
    L = accumarray(ic, triid, [], @(x){x});
    AdjT = sparse(Mt,Mt);
    for i=1:numel(L)
        T = L{i};
        if numel(T) >= 2
            for a=1:numel(T)
                for b=a+1:numel(T)
                    AdjT(T(a), T(b)) = 1;
                    AdjT(T(b), T(a)) = 1;
                end
            end
        end
    end
    AdjT = spones(AdjT);
end

function sigma = local_labels_to_sigma(labels, groups_order, params, clip)
% 0 -> SoftTissue, j>0 -> groups_order{j}
    N = numel(labels); sigma = zeros(N,1);
    cond = params.cond;
    if ~isfield(cond,'SoftTissue'), cond.SoftTissue = 0.35; end
    for j=1:numel(groups_order)
        gname = groups_order{j};
        if ~isfield(cond, gname), cond.(gname) = cond.SoftTissue; end
    end
    for i=1:N
        if labels(i)==0
            sigma(i) = cond.SoftTissue;
        else
            sigma(i) = cond.(groups_order{labels(i)});
        end
    end
    if ~isempty(clip) && all(isfinite(clip))
        sigma = min(max(sigma, clip(1)), clip(2));
    end
end

function L = local_levels_from_cond(condMap, groups_order)
% struct de niveaux (σ par groupe) pour logging
    L = struct();
    if isfield(condMap,'SoftTissue')
        L.SoftTissue = condMap.SoftTissue;
    else
        L.SoftTissue = 0.35;
    end
    for j=1:numel(groups_order)
        gj = groups_order{j};
        if isfield(condMap, gj)
            L.(gj) = condMap.(gj);
        end
    end
end

%% ===================== HELPERS GÉNÉRIQUES ===============================

function v = local_def(S, f, d)
    if ~isfield(S,f) || isempty(S.(f)), v = d; else, v = S.(f); end
end

function xy2 = local_simplify_ring(xy, step_mm)
% Ré-échantillonne une boucle à pas ~ step_mm (lissage + réduction de points).
    if isempty(xy), xy2 = xy; return; end
    if ~isequal(xy(1,:), xy(end,:)), xy = [xy; xy(1,:)]; end
    s = [0; cumsum(vecnorm(diff(xy),2,2))];
    L = s(end);
    if L <= 0, xy2 = xy(1:min(3,end),:); return; end
    M  = max(12, ceil(L / max(step_mm, 1e-6)));
    su = linspace(0, L, M).';
    xy2 = interp1(s, xy, su, 'linear');
end

function xy2 = local_clean_ring(xy, tol)
% Nettoie doublons consécutifs et boucle fermée propre.
    if nargin<2, tol = 0; end
    xy = xy(~any(isnan(xy),2), :);
    if ~isequal(xy(1,:), xy(end,:)), xy = [xy; xy(1,:)]; end
    d = hypot(diff(xy(:,1)), diff(xy(:,2)));
    keep = [true; d > max(tol, eps)];
    xy = xy(keep, :);
    if size(xy,1) >= 2 && isequal(xy(1,:), xy(end,:)), xy(end,:) = []; end
    xy2 = xy;
end
