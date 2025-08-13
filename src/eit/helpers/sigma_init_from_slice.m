function sigma0 = sigma_init_from_slice(results_src_mat, g_tgt, domain_tgt, soft_value)
% Projette une carte sigma (recon ou nodale) de la slice source -> noeuds de la slice cible.
S = load(results_src_mat);

% 1) carte nodale source
if isfield(S,'sigma_rec') && numel(S.sigma_rec)==size(S.g,1)
    sig_src = S.sigma_rec;
elseif isfield(S,'sigma_used') && numel(S.sigma_used)==size(S.g,1)
    sig_src = S.sigma_used;
elseif isfield(S,'sigma_tri') && numel(S.sigma_tri)==size(S.H,1)
    Nsrc = size(S.g,1);
    sig_sum = accumarray(S.H(:), repmat(S.sigma_tri,3,1), [Nsrc,1], @sum, 0);
    cnt     = accumarray(S.H(:), 1,                         [Nsrc,1], @sum, 0);
    sig_src = sig_sum ./ max(cnt,1);
    bad = ~isfinite(sig_src) | cnt==0;
    if any(bad)
        Ffill = scatteredInterpolant(S.g(~bad,1), S.g(~bad,2), sig_src(~bad), 'nearest', 'nearest');
        sig_src(bad) = Ffill(S.g(bad,1), S.g(bad,2));
    end
else
    error('Pas de carte sigma compatible dans %s', results_src_mat);
end

% 2) interpolation vers noeuds cibles
F = scatteredInterpolant(S.g(:,1), S.g(:,2), sig_src, 'natural', 'nearest');
sigma0 = F(g_tgt(:,1), g_tgt(:,2));

% 3) nettoyage / domaine
inside = isinterior(domain_tgt, g_tgt(:,1), g_tgt(:,2));
sigma0(~inside) = soft_value;
sigma0(~isfinite(sigma0)) = soft_value;
end
