function P2 = rmholes(P1)
% Retire les trous parasites d'un polyshape
    if P1.NumHoles == 0, P2 = P1; return; end
    regs = regions(P1); solids = regs([regs.IsSolid]);
    if isempty(solids), P2 = polyshape(); return; end
    P2 = solids(1);
    for k = 2:numel(solids), P2 = union(P2, solids(k)); end
end