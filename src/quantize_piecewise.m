function sigma_q = quantize_piecewise(sigma_nod, g, H, qopt)
% Quantifie sigma_nod (nodal) sur quelques niveaux constants.
% qopt.centers=[] => auto (k-means) avec qopt.k classes.
% qopt.clean_iters : nb cycles (lisse léger -> snap).

if ~isfield(qopt,'clean_iters'), qopt.clean_iters = 1; end
if ~isfield(qopt,'k'),           qopt.k = 4;          end
if ~isfield(qopt,'centers'),     qopt.centers = [];   end

Nn = size(g,1);

% valeurs par triangle (plus robuste pour clusteriser)
triVals = (sigma_nod(H(:,1)) + sigma_nod(H(:,2)) + sigma_nod(H(:,3)))/3;

% centres
if ~isempty(qopt.centers)
    centers = sort(qopt.centers(:));
else
    K = max(2, qopt.k);
    [~,C] = kmeans(triVals, K, 'Replicates',3, 'Start','plus');
    centers = sort(C(:));
end

% snap triangles -> centres
[~,ix] = min(abs(triVals - centers.'), [], 2);
triQ   = centers(ix);

% retour en nodal par moyenne des triangles incidents
sigma_q = tri2node_avg(H, triQ, Nn);

% nettoyage optionnel
for it=1:max(0,qopt.clean_iters)
    % lissage très léger
    L = mesh_laplacian(Nn, H);
    sigma_q = sigma_q - 0.2*(L*sigma_q);
    % re-snap nodal
    [~,ix] = min(abs(sigma_q - centers.'), [], 2);
    sigma_q = centers(ix);
end
end

% --- petits helpers ---
function nod = tri2node_avg(H, triVal, Nn)
acc = accumarray([H(:,1); H(:,2); H(:,3)], [triVal; triVal; triVal], [Nn,1], @mean, NaN);
if any(isnan(acc)), acc(isnan(acc)) = nanmean(acc); end
nod = acc;
end

function L = mesh_laplacian(Nn, H)
I = [H(:,1); H(:,1); H(:,2); H(:,2); H(:,3); H(:,3)];
J = [H(:,2); H(:,3); H(:,1); H(:,3); H(:,1); H(:,2)];
W = ones(numel(I),1);
A = sparse(I,J,W,Nn,Nn); A = A + A.'; A = A>0;
D = spdiags(sum(A,2),0,Nn,Nn);
L = D - A;
end
