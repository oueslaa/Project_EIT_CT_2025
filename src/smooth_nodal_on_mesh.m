function s = smooth_nodal_on_mesh(H, s, nIter, alpha, triGroup)
%SMOOTH_NODAL_ON_MESH  Petit lissage nodal par moyenne des voisins.
%   s = smooth_nodal_on_mesh(H, s, nIter, alpha)
%   s = smooth_nodal_on_mesh(H, s, nIter, alpha, triGroup)
%
% Entrées
%   H        : [Mt x 3] connectivité triangles (indices 1..N)
%   s        : [N x 1]  champ nodal à lisser (S/m)
%   nIter    : nb d'itérations (défaut 1)
%   alpha    : pas de lissage 0..1 (défaut 0.2)
%   triGroup : (optionnel) [Mt x 1] étiquette d'organe par triangle.
%              Si fourni, on NE lisse PAS à travers les frontières d'organes
%              (on déduit l'étiquette majoritaire de chaque nœud).
%
% Sortie
%   s        : champ nodal lissé (même taille que l'entrée)

    if nargin < 3 || isempty(nIter), nIter = 1; end
    if nargin < 4 || isempty(alpha), alpha = 0.2; end

    s = s(:);
    N = max(H(:));
    if numel(s) ~= N
        error('smooth_nodal_on_mesh: s doit être nodal (%d valeurs), reçu %d.', N, numel(s));
    end

    % --- construire la matrice d'adjacence nodale (bords non orientés)
    E = [H(:,[1 2]); H(:,[2 3]); H(:,[3 1])];
    A = sparse(E(:,1), E(:,2), 1, N, N);
    A = max(A, A.');                  % symétriser (u--v)

    % --- option: empêcher le lissage à travers les frontières d'organes
    if nargin >= 5 && ~isempty(triGroup)
        % étiquette majoritaire par nœud
        adjT = cell(N,1);
        for t = 1:size(H,1)
            v = H(t,:);
            adjT{v(1)}(end+1) = t; %#ok<AGROW>
            adjT{v(2)}(end+1) = t; %#ok<AGROW>
            adjT{v(3)}(end+1) = t; %#ok<AGROW>
        end
        nodeOrg = zeros(N,1);
        for i = 1:N
            tg = triGroup(adjT{i});
            nodeOrg(i) = mode(tg);
        end
        % annuler les arêtes qui relient des nœuds d'organes différents
        [ii, jj] = find(A);
        diffOrg = nodeOrg(ii) ~= nodeOrg(jj);
        A = A - sparse(ii(diffOrg), jj(diffOrg), 1, N, N);
        A = max(A, A.'); % resymétriser (au cas où)
    end

    deg = sum(A>0, 2);                % degré (nb voisins) par nœud
    nz  = deg > 0;

    for k = 1:nIter
        m = s;                         % par défaut, conserve la valeur si isolé
        m(nz) = (A(nz,:)*s) ./ deg(nz);
        s = (1-alpha)*s + alpha*m;
    end
end
