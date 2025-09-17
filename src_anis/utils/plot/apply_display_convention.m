function apply_display_convention(ax, mode)
% APPLY_DISPLAY_CONVENTION  Applique la convention d'affichage CT.
%
% But
% ----
% Harmoniser l'orientation des axes pour correspondre aux habitudes CT/DICOM.
%
% Entrées
% -------
% ax   : (optionnel) handle d'axes (par défaut: gca)
% mode : 'radiological' (gauche affichée à droite) ou 'neurological'
%        Défaut: 'radiological'
%
% Effets
% ------
% - Y vers le haut (coordonnées mm usuelles)
% - X inversé si 'radiological'
%
% Exemple
% -------
% imagesc(X,Y,I); apply_display_convention(gca,'radiological');
%
    if nargin==0 || isempty(ax), ax = gca; end
    if nargin<2 || isempty(mode), mode = 'radiological'; end
    set(ax,'YDir','normal');              % mm plots: Y up everywhere
    if strcmpi(mode,'radiological')
        set(ax,'XDir','reverse');         % match DICOM/CT display
    else
        set(ax,'XDir','normal');
    end
end
