function sm = fftSmooth(pts, K, step)
% FFTSMOOTH  Lisse un contour fermé par filtrage de Fourier.
%
% But
% ----
% Réduire le bruit / les petites irrégularités d'un contour 2D (Nx2, en mm
% ou pixels) en ne gardant que les K+1 premières fréquences de Fourier
% (symétriquement autour de la DC). Optionnellement, sous-échantillonne
% ensuite un point sur "step".
%
% Entrées
% -------
% pts  : [N x 2] points d'un contour fermé OU quasi-fermé (la fonction
%        referme si nécessaire). Les points doivent être ordonnés le long
%        du contour.
% K    : nombre de fréquences hautes à couper (K élevé = contour plus lisse).
% step : (optionnel) entier >=1. Si >1, prend 1 point sur "step" après IFFT.
%
% Sortie
% ------
% sm   : [M x 2] contour lissé (fermé : le 1er point = le dernier).
%
% Notes & pièges
% --------------
% - Si N<4, si pts contient des NaN/Inf, ou si tous les segments sont ~0,
%   la fonction renvoie le contour original (rien à lisser).
% - K est borné à floor((N-1)/2) pour respecter la symétrie hermitienne.
% - "unique(...,'stable')" évite les doublons exacts après l'IFFT.
% - La fermeture (duplication du 1er point en fin) est garantie.
%
% Exemple
% -------
% sm = fftSmooth(contour_mm, 50, 5);
%
    if size(pts,1)<4 || any(~isfinite(pts(:))), sm = pts; return; end
    d = hypot(diff(pts(:,1)), diff(pts(:,2))); if all(d<eps), sm = pts; return; end
    s = [0; cumsum(d)]/sum(d); N = numel(s);
    su = linspace(0,1,N)';
    Xu = interp1(s,pts(:,1),su,'linear',pts(1,1));
    Yu = interp1(s,pts(:,2),su,'linear',pts(1,2));
    FX = fft(Xu); FY = fft(Yu);
    K = min(K,floor((N-1)/2)); m = false(N,1); m([1:K+1, N-K+1:N])=true;
    FX(~m)=0; FY(~m)=0;
    pts2 = unique([real(ifft(FX)),real(ifft(FY))],'rows','stable');
    if nargin>2 && step>1, pts2=pts2(1:step:end,:); end
    if size(pts2,1)>=3 && ~isequal(pts2(1,:),pts2(end,:)), pts2(end+1,:)=pts2(1,:); end
    sm = pts2;
end
