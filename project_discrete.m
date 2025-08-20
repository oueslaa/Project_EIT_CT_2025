function sig_out = project_discrete(H, g, sig_in, levels, min_patch_tri)
% PROJECT_DISCRETE  Snap discret + nettoyage d'îlots sur un champ nodal.
% Inputs
%   H : (nTri x 3) connectivité
%   g : (nNode x 2/3) sommets (non utilisé ici mais gardé pour signature)
%   sig_in : (nNode x 1) champ nodal
%   levels : (k x 1) niveaux autorisés (S/m)
%   min_patch_tri : taille mini d'une composante (en triangles) à conserver
% Output
%   sig_out : (nNode x 1) champ nodal discretisé (valeurs ∈ levels)

if nargin<5 || isempty(min_patch_tri), min_patch_tri = 0; end
if isempty(levels) || numel(levels)==1
    sig_out = sig_in(:); return;
end
levels = levels(:);
nTri   = size(H,1);
nNode  = size(g,1);

% --- 1) Snap par triangle (sur moyenne nodale tri)
tri_val = (sig_in(H(:,1)) + sig_in(H(:,2)) + sig_in(H(:,3)))/3;
[~, tri_lab] = min(abs(tri_val - levels.'), [], 2);    % index de level par triangle
tri_discrete = levels(tri_lab);

% --- 2) Adjacence des triangles (partage d'arête)
A = tri_adjacency(H, nTri);

% --- 3) Suppression des petits îlots par label
if min_patch_tri > 0
    tri_lab = remove_small_patches(A, tri_lab, min_patch_tri, levels, tri_val);
    tri_discrete = levels(tri_lab);
end

% --- 4) Retour nodal : mode des labels voisins (vraiment discret)
node_ids = H(:);                       % 3*nTri x 1
lab_rep  = repelem(tri_lab, 3);        % 3*nTri x 1
nLev     = numel(levels);
counts   = accumarray([node_ids, lab_rep], 1, [nNode, nLev], @sum, 0);
[~, L]   = max(counts, [], 2);
sig_out  = levels(L);

end

% ================== helpers locales ===================

function A = tri_adjacency(H, nTri)
% Renvoie A (nTri x nTri) sparse, A(i,j)=1 si triangles i et j partagent une arête
E1 = sort(H(:,[1 2]), 2);
E2 = sort(H(:,[2 3]), 2);
E3 = sort(H(:,[1 3]), 2);
edges = [E1;E2;E3];
tidx  = [(1:nTri)'; (1:nTri)'; (1:nTri)'];

[edges_u, ~, ic] = unique(edges, 'rows');
% chaque arête unique apparaît 1 (bord) ou 2 fois (intérieur)
occ = accumarray(ic, 1);
pair_mask = occ(ic) == 2;                  % arêtes partagées
edges_sh  = edges(pair_mask, :);
tidx_sh   = tidx(pair_mask);

% retrouver les paires (t1,t2) pour chaque arête partagée
[~,~,ic2] = unique(edges_sh, 'rows');
g = accumarray(ic2, (1:numel(ic2))', [], @(v){v});  % indices des lignes par arête
I = []; J = [];
for k = 1:numel(g)
    idx = g{k};
    if numel(idx)==2
        t1 = tidx_sh(idx(1));  t2 = tidx_sh(idx(2));
        I = [I; t1; t2]; %#ok<AGROW>
        J = [J; t2; t1]; %#ok<AGROW>
    end
end
A = sparse(I, J, 1, nTri, nTri);
end

function tri_lab = remove_small_patches(A, tri_lab, min_sz, levels, tri_val)
% Pour chaque label, trouve les CC; les CC < min_sz sont fusionnées
nTri = size(A,1);
labels = unique(tri_lab(:).');
for L = labels
    mask = (tri_lab == L);
    if ~any(mask), continue; end

    % BFS/CC sur le sous-graphe A(mask,mask)
    idx_map = find(mask);
    subA = A(mask, mask);
    visited = false(nnz(mask),1);

    for s = 1:numel(idx_map)
        if visited(s), continue; end
        % BFS
        q = s; visited(s) = true; comp = s;
        while ~isempty(q)
            u = q(1); q(1) = [];
            neigh = find(subA(u,:));
            unvis = neigh(~visited(neigh));
            visited(unvis) = true;
            q = [q, unvis]; %#ok<AGROW>
            comp = [comp, unvis]; %#ok<AGROW>
        end
        comp_idx = idx_map(comp);          % indices globaux des triangles
        if numel(comp_idx) < min_sz
            % trouver les voisins hors composante
            N = find(any(A(comp_idx,:),1));            % tous voisins
            N = setdiff(N, comp_idx);
            if isempty(N)
                % îlot totalement isolé (théoriquement rare) -> on laisse
                continue;
            end
            neigh_lab = tri_lab(N);
            % label majoritaire autour; tie-break: niveau le + proche de la moyenne locale
            newL = mode(neigh_lab);
            if isnan(newL)
                m = mean(tri_val(comp_idx));
                [~, newL] = min(abs(levels - m));
            end
            tri_lab(comp_idx) = newL;
        end
    end
end
end
