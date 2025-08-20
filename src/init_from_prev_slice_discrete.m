function [sigma0_nod, info] = init_from_prev_slice_discrete(meshPrevFile, packPrevFile, meshThisFile, packThisFile, opts)
% INIT_FROM_PREV_SLICE_DISCRETE
%   Initialisation piecewise-constant depuis la slice "prev" vers la slice "this".
%   Retour: sigma0_nod (Nnodes_this x 1) STRICTEMENT quantifiée aux niveaux du forward.

if nargin < 5, opts = struct; end
clip   = getfielddef(opts, 'clip', []);
doWarp = getfielddef(opts, 'doWarp', true);

% --- chargement
MP = load(meshPrevFile);   % peut contenir 'contour'
PP = load(packPrevFile);   % g,H,sigma_tri,params.cond,...
MT = load(meshThisFile);   % g,H,contour
if nargin >= 4 && ~isempty(packThisFile) && isfile(packThisFile)
    PT = load(packThisFile); %#ok<NASGU>
end

gP = PP.g;  HP = PP.H;
gT = MT.g;  HT = MT.H;

% --- niveaux de conductivité (priorité params.cond)
levels = build_levels_from_prev(PP, getfielddef(opts,'init',struct()).quant);

% --- étiquettes des triangles de la prev (par proximité aux niveaux)
sigp_tri = double(PP.sigma_tri(:));
L        = double(levels(:).');
[~, lab_prev] = min(abs(sigp_tri - L), [], 2);   % lab_prev ∈ {1..K}

% --- centroïdes
cP  = (gP(HP(:,1),:) + gP(HP(:,2),:) + gP(HP(:,3),:)) / 3;
cT  = (gT(HT(:,1),:) + gT(HT(:,2),:) + gT(HT(:,3),:)) / 3;
cPw = cP; warp_used = false;

% --- warp contours si dispo (similarité)
if doWarp && isfield(MP,'contour') && isfield(MT,'contour') ...
          && ~isempty(MP.contour) && ~isempty(MT.contour)
    try
        N = min(size(MP.contour,1), size(MT.contour,1)); N = max(N, 50);
        C1s = resample_polyline(MP.contour,N);
        C2s = resample_polyline(MT.contour,N);
        tform = fitgeotrans(C1s, C2s, 'similarity');
        cPw  = transformPointsForward(tform, cP);
        warp_used = true;
    catch
        cPw = cP;
        warp_used = false;
    end
end

% --- transfert: pour chaque tri THIS, prendre le label du tri PREV le + proche
idx_prev_for_this = nearest_idx(cPw, cT);                % (nTriThis x 1)
lab_this = lab_prev(idx_prev_for_this);                  % (nTriThis x 1)

% --- sécurité dimension (au cas où)
nTriThis = size(HT,1);
if length(lab_this) < nTriThis
    lab_this(end+1:nTriThis,1) = mode(lab_this);
elseif length(lab_this) > nTriThis
    lab_this = lab_this(1:nTriThis);
end

% --- nettoyage (vote majoritaire)
Atri = tri_adjacency(HT);
clean_iters = getfielddef(getfielddef(opts,'init',struct()),'quant',struct());
clean_iters = getfielddef(clean_iters,'clean_iters',2);
for it = 1:max(0, clean_iters)
    lab_this = majority_step(lab_this, Atri);
end

% --- projection nodale par MODE (pas de moyenne)
Lnod = labels_to_nodes(lab_this, HT, size(gT,1));   % labels 1..K
sigma0_nod = levels( max(1, min(numel(levels), Lnod)) );   % snap exact via labels

% --- clip sécurité
if ~isempty(clip)
    sigma0_nod = min(max(sigma0_nod, clip(1)), clip(2));
end

% --- SNAP FINAL anti-valeurs intermédiaires (au cas où un traitement amont les a introduites)
sigma0_nod = snap_to_levels(sigma0_nod, levels);

% --- sortie
info = struct('levels',levels,'lab_tri',lab_this,'lab_nod',Lnod,'warp_used',warp_used);
end

% ==================== helpers ====================

function levels = build_levels_from_prev(PP, quant)
% Niveaux depuis params.cond si possible, sinon uniques de sigma_tri, sinon k-means/quantiles
if isfield(PP,'params') && isfield(PP.params,'cond') && isstruct(PP.params.cond)
    c = PP.params.cond;
    vals = [];
    fn = fieldnames(c);
    for i=1:numel(fn)
        vals(end+1) = c.(fn{i}); %#ok<AGROW>
    end
    vals = unique(vals);            % éviter duplicats (ex: Lung=Trachea)
    if ~isempty(vals), levels = sort(vals(:).'); return; end
end

levels = [];
if isfield(PP,'sigma_tri') && ~isempty(PP.sigma_tri)
    u = unique(PP.sigma_tri(:));
    if ~isempty(u), levels = sort(u(:).'); return; end
end

k = 4;
if isstruct(quant) && isfield(quant,'k') && ~isempty(quant.k), k = quant.k; end
st = PP.sigma_tri(:);
if exist('kmeans','file') == 2
    [~,C] = kmeans(st, k, 'Replicates',3);
    levels = sort(C(:).');
else
    qs = linspace(0,1,k+2); qs = qs(2:end-1);
    levels = sort(quantile(st, qs));
end
end

function v = getfielddef(S,f,d)
if isstruct(S) && isfield(S,f) && ~isempty(S.(f))
    v = S.(f);
else
    v = d;
end
end

function C = resample_polyline(P, N)
t = [0; cumsum(sqrt(sum(diff(P,1,1).^2,2)))];
if t(end) == 0
    C = repmat(P(1,:), N, 1); return;
end
t  = t / t(end);
tq = linspace(0,1,N).';
C  = [interp1(t, P(:,1), tq, 'linear'), interp1(t, P(:,2), tq, 'linear')];
end

function idx = nearest_idx(srcPts, dstPts)
if exist('knnsearch','file') == 2
    idx = knnsearch(srcPts, dstPts);
elseif exist('dsearchn','file') == 2
    idx = dsearchn(srcPts, dstPts);
else
    M = size(dstPts,1); N = size(srcPts,1);
    idx = zeros(M,1);
    for i=1:M
        d = sqrt(sum((srcPts - dstPts(i,:)).^2,2));
        [~,idx(i)] = min(d);
    end
end
end

function A = tri_adjacency(H)
E  = [H(:,[1 2]); H(:,[2 3]); H(:,[3 1])];
E  = sort(E,2);
[ue,~,ic] = unique(E, 'rows');
m = size(H,1);
A = sparse([],[],[], m, m, 6*m);
for e = 1:size(ue,1)
    tri = find(ic == e);
    if numel(tri) == 2
        A(tri(1), tri(2)) = 1;
        A(tri(2), tri(1)) = 1;
    end
end
A = spones(A);
end

function lab = majority_step(lab, A)
% Remplace chaque label par le mode de ses voisins (borne sécurité taille)
m = size(A,1); nLab = numel(lab);
m = min(m, nLab);
lab2 = lab;
for i = 1:m
    nb = find(A(i,:));
    nb = nb(nb <= nLab);
    if isempty(nb), continue; end
    L = lab(nb);
    lab2(i) = mode([lab(i); L(:)]);
end
lab = lab2;
end

function Lnod = labels_to_nodes(lab_tri, H, nNodes)
% Assigne à chaque nœud le MODE des labels des triangles incidents
Lnod = ones(nNodes,1) * mode(lab_tri); % fallback
cellT = accumarray(H(:), repmat((1:size(H,1)).',3,1), [nNodes,1], @(x){x});
for i = 1:nNodes
    tlist = cellT{i};
    if ~isempty(tlist)
        Lnod(i) = mode(lab_tri(tlist));
    end
end
end

function y = snap_to_levels(x, levels)
% Mappe chaque valeur de x au niveau le plus proche (discrétisation dure)
x = x(:);
L = levels(:).';
[~,idx] = min(abs(x - L), [], 2);
y = L(idx);
end
