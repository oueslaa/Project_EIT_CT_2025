classdef EITDiagnostics
% EITDiagnostics - Collecte & affichage des diagnostics (math + visuel)
% - SANS modifier SolverGN.m
% - Utilise le callback invSolver.plotter.plot(p) à chaque itération

methods(Static)

function ctx = new_context(args)
% args = struct with fields:
%   .stage_tag (char)    : 'A' / 'B' / 'custom'
%   .outdir   (char)     : dossier d'export (créé si absent)
%   .H, .g                mesh
%   .groupIdx (Mt x 1)   : mapping tri -> param
%   .y (m x 1)           : données ALIGNEES (après (a,b))
%   .L (m x m sparse)    : poids bruit
%   .lb, .ub (scalars)   : bornes
%   .p0 (K x 1)          : prior
%   .W  (K x K)          : régularisation (peut être speye)
%   .lambda (scalar)     : poids Tik
%   .solver              : EITFEM (pour SolveForwardVec)
%   .sigma_init_nodes(Nx1): pour delta
%   .save_iter_images (bool)
%   .save_metrics_csv (bool)
%   .realign_each_iter (bool) [optionnel, par défaut false]
%   .extra_plot_fcn (function handle) [optionnel] pour dessiner aussi les maps iteratives

    arguments
        args struct
    end
    ctx = args;
    if ~isfield(ctx,'realign_each_iter'), ctx.realign_each_iter = false; end
    if ~isfield(ctx,'extra_plot_fcn'),    ctx.extra_plot_fcn    = [];    end
    if ~isempty(ctx.outdir) && ~exist(ctx.outdir,'dir'), mkdir(ctx.outdir); end

    % historiques
    ctx.iter        = 0;
    ctx.p_prev      = [];
    ctx.p_hist      = {};           % cell per iter
    ctx.misfit_norm = [];           % ||L r||2
    ctx.reg_term    = [];           % lambda * ||W (p-p0)||^2
    ctx.obj_val     = [];           % 0.5*misfit^2 + 0.5*reg
    ctx.step_norm   = [];           % ||p_k - p_{k-1}||
    ctx.a_hist      = [];           % si realign_each_iter = true
    ctx.b_hist      = [];

    % stockage rapide des dernières cartes pour le dashboard
    ctx.last_sig_nodes = [];
    ctx.last_delta     = [];
end

function iter_and_plot(p, ctx)
% Wrapper : log + (optionnel) affiche la carte "propre" de l'itération
    EITDiagnostics.iter_callback(p, ctx);
    if ~isempty(ctx.extra_plot_fcn)
        try, ctx.extra_plot_fcn(p); end %#ok<TRYNC>
    end
end

function iter_callback(p, ctx)
% Appelé à CHAQUE itération par GN (via plotter.plot)
% - calcule U(p), résidu, termes misfit/reg, norme du pas, (a,b) optionnels
% - exporte une image de la carte si demandé
    ctx.iter = ctx.iter + 1;
    H = ctx.H; g = ctx.g; groupIdx = ctx.groupIdx;
    y = ctx.y; L = ctx.L; lb = ctx.lb; ub = ctx.ub;
    p0 = ctx.p0; W = ctx.W; lambda = ctx.lambda;

    % Param → sigma_tri → tensions simulées
    sig_tri = p(groupIdx);
    U = ctx.solver.SolveForwardVec(sig_tri); U = U(:);

    % Ré-alignement optionnel (robuste), sinon y est déjà aligné
    a = 1; b = 0;
    if ctx.realign_each_iter
        A = [U, ones(numel(U),1)];
        ab = A \ y; a = ab(1); b = ab(2);
        y_eff = (y - b) / max(a, eps);
    else
        y_eff = y;
    end

    r = y_eff - U;
    misfit = norm(L*r);
    diffp  = p - p0;
    regv   = lambda * sum((W*diffp).^2);
    obj    = 0.5*misfit^2 + 0.5*regv;

    steplen = NaN;
    if ~isempty(ctx.p_prev)
        steplen = norm(p - ctx.p_prev);
    end
    ctx.p_prev = p;

    % Nodal pour visuels
    N = size(g,1);
    sig_nodes = EITDiagnostics.tri2node_avg(H, sig_tri, N);
    sig_nodes = min(max(sig_nodes, lb), ub);
    delta = sig_nodes - ctx.sigma_init_nodes;

    % Historique
    ctx.p_hist{end+1}      = p;
    ctx.misfit_norm(end+1) = misfit;
    ctx.reg_term(end+1)    = regv;
    ctx.obj_val(end+1)     = obj;
    ctx.step_norm(end+1)   = steplen;
    ctx.a_hist(end+1)      = a;
    ctx.b_hist(end+1)      = b;
    ctx.last_sig_nodes     = sig_nodes;
    ctx.last_delta         = delta;

    % Export image d'itération (léger)
    if isfield(ctx,'save_iter_images') && ctx.save_iter_images && ~isempty(ctx.outdir)
        f = figure('Visible','off','Color','w');
        t = tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');
        ax1 = nexttile(t,1); ax2 = nexttile(t,2);
        EITDiagnostics.draw_field(ax1, H, g, sig_nodes); axis(ax1,'equal'); axis(ax1,'tight');
        colormap(ax1, turbo); caxis(ax1,[lb ub]); colorbar(ax1); title(ax1, sprintf('\\sigma (iter %d)', ctx.iter),'Interpreter','tex');
        EITDiagnostics.draw_field(ax2, H, g, delta); axis(ax2,'equal'); axis(ax2,'tight');
        colormap(ax2, parula); Lc = 0.5*(ub-lb); caxis(ax2,[-Lc Lc]); colorbar(ax2); title(ax2,'\Delta\sigma','Interpreter','tex');
        exportgraphics(f, fullfile(ctx.outdir, sprintf('iter_%03d.png', ctx.iter)), 'Resolution', 160, 'BackgroundColor','white');
        close(f);
    end
end

function finish(ctx)
% Construit un DASHBOARD + export CSV des métriques
    it = (1:numel(ctx.obj_val))';
    T = table(it, ctx.misfit_norm(:), ctx.reg_term(:), ctx.obj_val(:), ctx.step_norm(:), ctx.a_hist(:), ctx.b_hist(:), ...
              'VariableNames', {'iter','misfit_norm','reg_term','obj','step_norm','a','b'});

    if isfield(ctx,'save_metrics_csv') && ctx.save_metrics_csv && ~isempty(ctx.outdir)
        writetable(T, fullfile(ctx.outdir,'metrics.csv'));
    end

    % --- Fallback : reconstruit la carte finale si besoin ---
    if isempty(ctx.last_sig_nodes)
        p_final = ctx.p_hist{end};
        sig_tri = p_final(ctx.groupIdx);
        N = size(ctx.g,1);
        sig_nodes = EITDiagnostics.tri2node_avg(ctx.H, sig_tri, N);
        sig_nodes = min(max(sig_nodes, ctx.lb), ctx.ub);
        ctx.last_sig_nodes = sig_nodes;
        ctx.last_delta     = sig_nodes - ctx.sigma_init_nodes;
    end

    % Dashboard (6 panneaux)
    f = figure('Color','w','Name',sprintf('Diagnostics Stage %s', ctx.stage_tag));
    t = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');

    % (1) sigma_init
    ax1 = nexttile(t,1);
    EITDiagnostics.draw_field(ax1, ctx.H, ctx.g, ctx.sigma_init_nodes);
    axis(ax1,'equal'); axis(ax1,'tight'); colormap(ax1,turbo); caxis(ax1,[ctx.lb ctx.ub]); colorbar(ax1);
    title(ax1,'\sigma_{init}','Interpreter','tex');

    % (2) sigma_final
    ax2 = nexttile(t,2);
    EITDiagnostics.draw_field(ax2, ctx.H, ctx.g, ctx.last_sig_nodes);
    axis(ax2,'equal'); axis(ax2,'tight'); colormap(ax2,turbo); caxis(ax2,[ctx.lb ctx.ub]); colorbar(ax2);
    title(ax2,'\sigma_{final}','Interpreter','tex');

    % (3) delta final
    ax3 = nexttile(t,3);
    EITDiagnostics.draw_field(ax3, ctx.H, ctx.g, ctx.last_delta);
    axis(ax3,'equal'); axis(ax3,'tight'); colormap(ax3,parula); Lc=0.5*(ctx.ub-ctx.lb); caxis(ax3,[-Lc Lc]); colorbar(ax3);
    title(ax3,'\Delta\sigma_{final}','Interpreter','tex');

    % (4) Convergence (handles explicites pour éviter le warning de legend)
    ax4 = nexttile(t,4);
    h1 = plot(ax4, it, ctx.obj_val, '-o'); hold(ax4,'on');
    h2 = plot(ax4, it, ctx.misfit_norm, '-s');
    yyaxis(ax4,'right');
    h3 = plot(ax4, it, ctx.step_norm, '-^');
    yyaxis(ax4,'left');
    grid(ax4,'on'); xlabel(ax4,'iter'); ylabel(ax4,'obj / ||Lr||');
    title(ax4,'Convergence');
    legend([h1 h2 h3], {'obj','||Lr||','||\Delta p||'}, 'Location','best');

    % (5) Scatter U vs y (final)
    p_final = ctx.p_hist{end};
    sig_tri = p_final(ctx.groupIdx);
    U = ctx.solver.SolveForwardVec(sig_tri); U = U(:);
    y_eff = ctx.y;
    if ctx.realign_each_iter
        A = [U, ones(numel(U),1)]; ab = A \ y_eff; a=ab(1); b=ab(2);
        y_eff = (y_eff - b)/max(a,eps);
    end
    ax5 = nexttile(t,5);
    plot(ax5, y_eff, U, '.'); grid(ax5,'on'); xlabel(ax5,'y (aligné)'); ylabel(ax5,'U(\sigma_{final})');
    % diagonale
    mn = min(min(y_eff), min(U)); mx = max(max(y_eff), max(U));
    hold(ax5,'on'); plot(ax5, [mn mx], [mn mx], '-'); hold(ax5,'off');
    title(ax5,'Scatter y vs U (final)');

    % (6) Histogramme résidu final
    r = y_eff - U;
    ax6 = nexttile(t,6);
    histogram(ax6, r, 'Normalization','pdf');
    grid(ax6,'on'); title(ax6,'Résidu final (pdf)'); xlabel(ax6,'r_i'); ylabel(ax6,'pdf');

    if ~isempty(ctx.outdir)
        exportgraphics(f, fullfile(ctx.outdir,'dashboard.png'), 'Resolution', 220, 'BackgroundColor','white');
    end
end

% ---- Petites utilitaires visuelles ----
function draw_field(ax, H, g, C)
    if nargin<1 || isempty(ax), ax=gca; end
    C=C(:); N=size(g,1); M=size(H,1); hold(ax,'on');
    if numel(C)==N
        patch('Faces',H,'Vertices',g,'FaceVertexCData',C,'FaceColor','interp', ...
              'EdgeColor','none','Parent',ax);
    elseif numel(C)==M
        patch('Faces',H,'Vertices',g,'FaceVertexCData',C,'FaceColor','flat', ...
              'EdgeColor','none','Parent',ax);
    else
        Cnod = EITDiagnostics.tri2node_avg(H, C, N);
        patch('Faces',H,'Vertices',g,'FaceVertexCData',Cnod,'FaceColor','interp', ...
              'EdgeColor','none','Parent',ax);
    end
end

function nod = tri2node_avg(H, tri, N)
    tri = tri(:); Mt=size(H,1);
    if numel(tri) ~= Mt
        if mod(numel(tri),Mt)==0
            tri = mean(reshape(tri,[],Mt).',2);
        else
            tri = tri(1:Mt);
        end
    end
    nod = zeros(N,1); cnt = zeros(N,1);
    for t=1:Mt
        v = H(t,:);
        nod(v) = nod(v) + tri(t);
        cnt(v) = cnt(v) + 1;
    end
    idx = cnt>0; nod(idx) = nod(idx)./cnt(idx);
end

end % methods
end
