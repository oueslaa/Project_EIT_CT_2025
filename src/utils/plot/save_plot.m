function paths = save_plot(fig, outfile_noext, opts)
% SAVE_PLOT  Export a figure to PNG (+ optional PDF, SVG, FIG).
%   paths = save_plot(fig, 'Outputs/s0011/slice_301/sigma_forward', opts)
%   opts fields (all optional):
%     .dpi  (default 220)   .pdf (true/false)   .svg (false)   .fig (false)
%     .mode (radiological/neurological) applied via apply_display_convention
%
% Returns struct with written paths.

    if nargin<3, opts = struct; end
    if ~isfield(opts,'dpi'),  opts.dpi = 220; end
    if ~isfield(opts,'pdf'),  opts.pdf = true; end
    if ~isfield(opts,'svg'),  opts.svg = false; end
    if ~isfield(opts,'fig'),  opts.fig = false; end
    if ~isfield(opts,'mode'), opts.mode = 'radiological'; end

    outfile_noext = char(outfile_noext);
    outdir = fileparts(outfile_noext);
    if ~exist(outdir,'dir'), mkdir(outdir); end

    ax = gca;
    % Orientation + toolbar
    try, apply_display_convention(ax, opts.mode); end %#ok<TRYNC>
    try
        if isprop(ax,'Toolbar') && ~isempty(ax.Toolbar) && isprop(ax.Toolbar,'Visible')
            ax.Toolbar.Visible = 'off';
        end
    end %#ok<TRYNC>

    % PNG
    png = [outfile_noext '.png'];
    exportgraphics(fig, png, 'Resolution', opts.dpi, 'BackgroundColor','white');
    paths.png = png;

    % Vector (PDF/SVG)
    if opts.pdf
        pdf = [outfile_noext '.pdf'];
        exportgraphics(fig, pdf, 'ContentType','vector', 'BackgroundColor','white');
        paths.pdf = pdf;
    end
    if opts.svg
        svg = [outfile_noext '.svg'];
        exportgraphics(fig, svg, 'ContentType','vector', 'BackgroundColor','white');
        paths.svg = svg;
    end

    % FIG (editable)
    if opts.fig
        figf = [outfile_noext '.fig'];
        savefig(fig, figf);
        paths.fig = figf;
    end
end
