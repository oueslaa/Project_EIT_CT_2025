function tri = node2tri_avg_local(H, nod)
nod = nod(:);
tri = (nod(H(:,1))+nod(H(:,2))+nod(H(:,3)))/3;
tri = tri(:);
end
