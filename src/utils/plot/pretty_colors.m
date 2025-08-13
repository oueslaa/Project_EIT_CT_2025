function cols = pretty_colors(n)
base = [0.2 0.6 1; 0.9 0.4 0.1; 1 0.8 0; 0.6 0 0.8; 0.2 0.8 0; 0.3 0.3 0.3];
if n <= size(base,1), cols = base(1:n,:); else, cols = lines(n); end
end