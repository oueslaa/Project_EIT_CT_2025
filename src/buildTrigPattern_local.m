function I = buildTrigPattern_local(Ne, Iamp)
I = zeros(Ne, Ne);
for k = 1:Ne
    ph = 2*pi*(k-1)/Ne;
    for n = 1:Ne
        I(n, k) = Iamp * sin(2*pi*(n-1)/Ne + ph);
    end
end
end
