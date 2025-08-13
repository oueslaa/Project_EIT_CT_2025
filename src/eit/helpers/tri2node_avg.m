function sigma_nod = tri2node_avg(H, sigma_tri, Nnodes)
% Moyenne des valeurs triangulaires vers les nœuds
    acc = accumarray(H(:), repmat(sigma_tri,3,1), [Nnodes,1], @sum, 0);
    cnt = accumarray(H(:), 1, [Nnodes,1], @sum, 0);
    sigma_nod = acc ./ max(cnt,1);
end