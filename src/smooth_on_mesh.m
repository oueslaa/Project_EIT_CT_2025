function s = smooth_on_mesh(s, H, lam, iters)
% s : vecteur nodal (longueur = nb de noeuds)
% H : connectivité triangles (nT x 3)
% lam in [0,1], iters >= 0

    % -- défauts
    if nargin < 4 || isempty(iters), iters = 1; end
    if nargin < 3 || isempty(lam),   lam   = 0.20; end

    % -- vecteur colonne + types sûrs
    s     = double(s(:));
    lam   = double(lam);
    iters = double(iters);

    % -- taille attendue
    N = max(H(:));             % nb de noeuds attendus par H
    if numel(s) ~= N
        % Si s est homogène mais mauvaise forme, essaye d'adapter prudemment
        % (mais en EIT normalement numel(s) == max(H(:))).
        error('smooth_on_mesh: size mismatch: numel(s)=%d, max(H)=%d.', numel(s), N);
    end

    % -- matrice d'adjacence (umbrella)
    i1 = H(:,1); i2 = H(:,2); i3 = H(:,3);
    A  = sparse([i1;i2;i3],[i2;i3;i1], 1, N, N);
    A  = A + A.';                                % symétrise
    A  = A - spdiags(diag(A),0,N,N);             % pas de boucle sur soi-même

    deg = full(sum(A,2)); 
    deg(deg==0) = 1;
    Dinv = spdiags(1./deg,0,N,N);
    W    = Dinv * A;                             % moyenne des voisins

    % -- itérations de lissage
    for t = 1:iters
        s = (1-lam)*s + lam*(W*s);
    end
end
