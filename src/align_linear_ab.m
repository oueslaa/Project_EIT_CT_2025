function [a,b,pre,post] = align_linear_ab(U, M)
% Trouve a,b tel que || (M-b)/a - U ||_2 soit minimale.
U = U(:); M = M(:);
Um = mean(U); Mm = mean(M);
Ud = U - Um; Md = M - Mm;
den = dot(Ud,Ud);
if den < eps
    a = 1; b = 0;
else
    a = dot(Ud, Md) / den;
    b = Mm - a*Um;
end
pre  = norm(M - U);
post = norm((M - b)/max(a,eps) - U);
end
