function apply_display_convention(ax, mode)
% mode: 'radiological' (left appears on the right) or 'neurological'
    if nargin==0 || isempty(ax), ax = gca; end
    if nargin<2 || isempty(mode), mode = 'radiological'; end
    set(ax,'YDir','normal');              % mm plots: Y up everywhere
    if strcmpi(mode,'radiological')
        set(ax,'XDir','reverse');         % match DICOM/CT display
    else
        set(ax,'XDir','normal');
    end
end
