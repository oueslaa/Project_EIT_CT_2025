function triVals = node2tri_avg(H, nodal)
% Moyenne nodale -> valeur par triangle
nodal = nodal(:);
triVals = (nodal(H(:,1)) + nodal(H(:,2)) + nodal(H(:,3)))/3;
end
