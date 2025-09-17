function P2 = rmholes(P1)
% RMHOLES  Supprime les trous "parasites" d'un polyshape.
%
% But
% ----
% Nettoyer un polyshape pour ne conserver que les régions solides (sans trous).
%
% Entrée
% ------
% P1 : polyshape quelconque (éventuellement multipolygone).
%
% Sortie
% ------
% P2 : polyshape constitué de l'union de toutes les régions solides de P1.
%
% Exemple
% -------
% body_clean = rmholes(body_raw);
%
    if P1.NumHoles == 0, P2 = P1; return; end
    regs = regions(P1); solids = regs([regs.IsSolid]);
    if isempty(solids), P2 = polyshape(); return; end
    P2 = solids(1);
    for k = 2:numel(solids), P2 = union(P2, solids(k)); end
end
