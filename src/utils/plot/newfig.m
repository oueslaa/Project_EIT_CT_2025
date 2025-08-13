function fig = newfig(name, visible)
% NEWFIG  Create a figure with sensible defaults for batch/export.
%   fig = newfig('Title', 'off');
    if nargin<2 || isempty(visible), visible = 'off'; end
    fig = figure('Name',name, 'NumberTitle','off', ...
                 'Color','w', 'Visible',visible);
    % Optional: unify look & feel (lightweight defaults)
    ax = axes('Parent',fig); %#ok<LAXES>
    set(ax, 'Box','on', 'Layer','top', 'FontName','Helvetica');
end
