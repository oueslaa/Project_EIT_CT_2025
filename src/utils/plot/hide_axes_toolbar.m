function hide_axes_toolbar(ax)
% Hides the axes toolbar if present (MATLAB R2019b+)
    if nargin==0 || isempty(ax), ax = gca; end
    try
        if isprop(ax,'Toolbar') && ~isempty(ax.Toolbar) && isprop(ax.Toolbar,'Visible')
            ax.Toolbar.Visible = 'off';
        end
    catch
        % older MATLAB: ignore silently
    end
end
