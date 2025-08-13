function P = makeSimplePlotterEIT(g, H)
% Plotter minimal compatible SolverGN (fields .plot, .plotUQ, .lineplot)
P.plot     = @(sigma) local_plot_sigma(g,H,sigma);
P.plotUQ   = @(w) [];
P.lineplot = @(sigma, uq, t) [];
end

function local_plot_sigma(g,H,sigma)
persistent fh ax
if isempty(fh) || ~isgraphics(fh)
    fh = figure('Name','GN: sigma (iteration)','Color','w');
    ax = axes('Parent',fh); hold(ax,'on'); axis(ax,'equal'); axis(ax,'tight');
    colormap(ax, jet);
else
    figure(fh); axes(ax); cla(ax);
end
C = sigma(:);
N = size(g,1); M = size(H,1);
if numel(C)==N
    patch('Faces',H,'Vertices',g,'FaceVertexCData',C,'FaceColor','interp','EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',0.3,'Parent',ax);
elseif numel(C)==M
    patch('Faces',H,'Vertices',g,'FaceVertexCData',C,'FaceColor','flat','EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',0.3,'Parent',ax);
else
    % fallback tri->noeud
    Cn = zeros(N,1); cnt=zeros(N,1);
    for t=1:size(H,1)
        v = H(t,:);
        val = C(min(t,end));
        Cn(v) = Cn(v)+val; cnt(v)=cnt(v)+1;
    end
    Cn(cnt>0)=Cn(cnt>0)./cnt(cnt>0);
    patch('Faces',H,'Vertices',g,'FaceVertexCData',Cn,'FaceColor','interp','EdgeColor',[0.8 0.8 0.8],'EdgeAlpha',0.3,'Parent',ax);
end
colorbar(ax); title(ax,'\sigma estimate (GN)','Interpreter','tex');
set(ax,'XTick',[],'YTick',[]); drawnow;
end
