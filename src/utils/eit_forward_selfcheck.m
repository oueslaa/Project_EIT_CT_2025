function stats = eit_forward_selfcheck(rootOut)
% stats = eit_forward_selfcheck('Outputs/s0011/slice_302')
S = load(fullfile(rootOut,'eit_pack.mat'));   % contains Imeas, Imeas_clean, params
Ne = S.params.Ne;
V  = reshape(S.Imeas_clean, Ne, []);          % propre
I  = EITSim.buildTrigPattern(Ne, S.params.I_amp);

% 1) somme des courants (chaque motif doit sommer à 0)
sumI = sum(I,1);
stats.sumI_maxabs = max(abs(sumI));

% 2) réciprocité sans inversion (symétrie de I^T V)
Smat = I.' * V;
stats.recip_error = norm(Smat - Smat.', 'fro') / norm(Smat, 'fro');

% 3) corrélation colonne k vs k+1 décalée d'une électrode
K = size(V,2);
c = zeros(1, K-1);
for k=1:K-1
    vk  = V(:,k);
    vkp = circshift(V(:,k+1),1);
    c(k) = (vk.'*vkp) / (norm(vk)*norm(vkp));
end
stats.shift_corr_mean = mean(c);

% 4) ordre de grandeur (mV)
stats.V_pk_mV = 1e3*max(abs(V(:)));

% print
fprintf('sum(I) max abs: %.2e | reciprocity: %.3e | shift corr: %.3f | Vpk ~ %.2f mV\n', ...
    stats.sumI_maxabs, stats.recip_error, stats.shift_corr_mean, stats.V_pk_mV);
end
