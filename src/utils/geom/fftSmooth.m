function sm = fftSmooth(pts, K, step)
% Lisse un contour fermé (Nx2) par filtrage de Fourier
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