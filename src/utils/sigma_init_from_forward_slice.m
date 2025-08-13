function sigma0 = sigma_init_from_forward_slice(res_mat_src, g_tgt, domain_tgt, soft_val)
% SIGMA_INIT_FROM_FORWARD_SLICE
%   Construit sigma_init sur les NOEUDS de la cible en prenant la carte
%   forward de la slice source (sigma_tri) et en la projetant par
%   localisation de points dans le maillage source.
%
% Inputs
%   - res_mat_src : Outputs/.../slice_XXX/eit_results.mat de la slice source
%   - g_tgt       : [Nt x 2] noeuds (mm) de la cible
%   - domain_tgt  : polyshape (mm) cible, pour sécuriser l'extérieur
%   - soft_val    : valeur SoftTissue (S/m) pour le remplissage
%
% Output
%   - sigma0      : [Nt x 1] carte nodale (S/m) pour sigma_init

S = load(res_mat_src);
assert(isfield(S,'g') && isfield(S,'H'), ...
    'Le fichier source doit contenir g et H.');
g_src = S.g; H_src = S.H;

% --- on veut PRIORITAIREMENT la carte forward triangulaire ---
if isfield(S,'sigma_tri') && numel(S.sigma_tri)==size(H_src,1)
    sig_tri = S.sigma_tri(:);
elseif isfield(S,'sigma_used') && numel(S.sigma_used)==size(H_src,1)
    % selon certains solveurs, "sigma_used" peut être triangulaire
    sig_tri = S.sigma_used(:);
else
    error('Aucune carte triangulaire trouvée (sigma_tri/sigma_used) dans %s.', res_mat_src);
end

% Triangulation de la slice source
TR = triangulation(H_src, g_src);

% Localisation des noeuds cibles dans les triangles source
tid = pointLocation(TR, g_tgt);   % indice de triangle ou NaN
sigma0 = soft_val * ones(size(g_tgt,1),1);
ok = ~isnan(tid);
sigma0(ok) = sig_tri(tid(ok));

% Pour les noeuds non localisés (NaN), projeter sur le triangle le plus proche
if any(~ok)
    % centres des triangles source
    Csrc = (g_src(H_src(:,1),:) + g_src(H_src(:,2),:) + g_src(H_src(:,3),:))/3;
    idx_near = knnsearch(Csrc, g_tgt(~ok,:), 'K', 1);
    sigma0(~ok) = sig_tri(idx_near);
end

% Sécurité : tout ce qui est en dehors du domaine cible -> SoftTissue
try
    in = isinterior(domain_tgt, g_tgt(:,1), g_tgt(:,2));
    sigma0(~in) = soft_val;
end
end
