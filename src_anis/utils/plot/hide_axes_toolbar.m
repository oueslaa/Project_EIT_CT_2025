function hide_axes_toolbar(ax)
% HIDE_AXES_TOOLBAR  Cache la mini-toolbar sur les axes (R2019b+).
%
% Entrée
% ------
% ax : (optionnel) axes cible (défaut: gca)
%
% Notes
% -----
% - Sans effet sur anciennes versions MATLAB (catch silencieux).
%
% Exemple
% -------
% hide_axes_toolbar(gca);
%
    if nargin==0 || isempty(ax), ax = gca; end
    try
        if isprop(ax,'Toolbar') && ~isempty(ax.Toolbar) && isprop(ax.Toolbar,'Visible')
            ax.Toolbar.Visible = 'off';
        end
    catch
        % older MATLAB: ignore silently
    end
end
