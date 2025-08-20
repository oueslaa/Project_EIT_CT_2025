function [sigma0_nod, info] = init_from_prev_slice(meshPrevFile, packPrevFile, meshThisFile, softDefault, opts)
% INIT_FROM_PREV_SLICE
%   Construit un sigma_init "par morceaux constants" pour la slice cible
%   en copiant la palette (valeurs constantes) de la slice précédente.
%
% OUT:
%   sigma0_nod  : conductivité nodale (constante par organe)
%   info        : struct debug (palette, labels, transform, etc.)

% --------- options (défauts + override) ----------
S.doWarp        = true;
S.clip          = [0, inf];
S.blend_homog   = 0.05;     % petit mélange vers l'homogène
S.palette_round = 1e-3;     % tolérance pour regrouper les niveaux constants
S.clean_iters   = 1;        % itérations de majorité sur triangles
if nargin >= 5 && ~isempty(opts) && isfield(opts,'init') && isstruct(opts.init)
    S = override_(S, opts.init);
end
if nargin >=5 && isfield(opts,'doWarp'), S.doWarp = opts.doWarp; end
if nargin >=5 && isfield(opts,'clip')   && ~isempty(opts.clip), S.clip = opts.clip; end

% --------- load fichiers ----------
MP = load(meshPrevFile);     % peut contenir MP.contour
PP = load(packPrevFile);     % gP, HP, sigma_tri (slice prev)
MT = load(meshThisFile);     % gT, HT (slice this)

% mailles
gP  = PP.g;   HP = PP.H;
gT  = MT.g;   HT = MT.H;
ntP = size(HP,1); ntT = size(HT,1); nT = size(gT,1);

% valeurs tri de la prev (palette source)
valsP = PP.sigma_tri(:);

% --------- (1) Warp des barycentres de la prev vers la géométrie courante
cP = (gP(HP(:,1),:) + gP(HP(:,2),:) + gP(HP(:,3),:)) / 3;  % (ntP x 2)
cT = (gT(HT(:,1),:) + gT(HT(:,2),:) + gT(HT(:,3),:)) / 3;  % (ntT x 2)
Tform = []; cPw = cP;

if S.doWarp && isfield(MP,'contour') && isfield(MT,'contour') ...
             && ~isempty(MP.contour) && ~isempty(MT.contour)
    Cprev = resample_polyline_(MP.contour, 200);
    Cthis = resample_polyline_(MT.contour, 200);
    try
        [~, Z, tr] = procrustes(Cthis, Cprev, 'Scaling', true, 'Reflection', false);
        % Z = Cprev transformé → approx Cthis
        % Applique la même transform aux points quelconques P:
        %   P*T.T*tr.b + tr.c
        Tform = tr;
        cPw   = cP*tr.T*tr.b + tr.c;
    catch
        % fallback = pas de warp
        Tform = [];
        cPw   = cP;
    end
end

% --------- (2) Palette des niveaux constants sur la prev
% Regroupement par arrondi (tolérance palette_round)
bins  = round(valsP / S.palette_round) * S.palette_round;
ubins = unique(bins);
K     = numel(ubins);              % nb de niveaux constants
labP  = zeros(ntP,1);
palette = zeros(K,1);
for k = 1:K
    in = (bins == ubins(k));
    labP(in)   = k;
    palette(k) = median(valsP(in)); % niveau "canonique"
end

% --------- (3) Transport des labels vers la nouvelle maille (NN)
idx = knnsearch(cPw, cT);        % pour chaque tri cible, le tri source le + proche
labT = labP(idx);                % labels de triangles sur la cible

% --------- (4) Nettoyage par majorité sur l'adjacence des triangles
if S.clean_iters > 0
    TR = triangulation(HT, gT);
    Nbr = TR.neighbors;          % (ntT x 3) indices tri voisins ou -1
    for it = 1:S.clean_iters
        labT_new = labT;
        for i = 1:ntT
            neigh = Nbr(i,:);  neigh = neigh(neigh>0);
            if ~isempty(neigh)
                m = mode([labT(i); labT(neigh)]);
                labT_new(i) = m;
            end
        end
        labT = labT_new;
    end
end

% --------- (5) Vote majoritaire "tri → nœuds" (pas de moyenne)
%   Pour chaque nœud, on prend le label majoritaire des triangles incidents
Smat = sparse(HT(:), repelem((1:ntT)',3,1), 1, nT, ntT);  % nœuds x triangles (incidence)
counts = zeros(nT, K);
for k = 1:K
    counts(:,k) = Smat * double(labT == k);
end
[~, labN] = max(counts, [], 2);   % label par nœud
sigma0_nod = palette(labN);

% --------- (6) Blend vers homogène + clip
if S.blend_homog > 0
    sigma0_nod = (1-S.blend_homog)*sigma0_nod + S.blend_homog*softDefault;
end
sigma0_nod = max(sigma0_nod, 1e-9);
if any(isfinite(S.clip))
    sigma0_nod = min(max(sigma0_nod, S.clip(1)), S.clip(2));
end

% --------- info debug
info = struct('palette', palette, 'K', K, 'lab_tri', labT, 'lab_node', labN, ...
              'tform', Tform, 'doWarp', S.doWarp);
end

% ===================== helpers =====================

function S = override_(S, U)
% remplace S.f par U.f si présent (récursif basique)
if ~isstruct(U), return; end
fn = fieldnames(U);
for i=1:numel(fn)
    f = fn{i};
    if isfield(S,f) && isstruct(S.(f)) && isstruct(U.(f))
        S.(f) = override_(S.(f), U.(f));
    else
        S.(f) = U.(f);
    end
end
end

function C = resample_polyline_(P, N)
% rééchantillonne une polyligne 2D en N pts
L = [0; cumsum(sqrt(sum(diff(P,1,1).^2,2)))];
if L(end) < eps, C = repmat(P(1,:), N, 1); return; end
t  = L / L(end);
tq = linspace(0,1,N)';
C  = [interp1(t, P(:,1), tq, 'linear'), ...
      interp1(t, P(:,2), tq, 'linear')];
end
