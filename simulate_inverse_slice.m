function rec = simulate_inverse_slice(patient_id, slice_this, slice_prev, opts)
% SIMULATE_INVERSE_SLICE — Reconstruction EIT (une passe GN + TV + L2-to-init)
% Simple, robuste, aucune sophistication inutile.
%
% Inputs:
%   patient_id, slice_this, slice_prev
%   opts : (tous facultatifs)
%       .clip [min max] (def [0.05 0.50])
%       .doWarp (def true)
%       .init.discrete (def true)
%       .init.calib_enable (def true)
%       .init.calib_clip [0.95 1.05]
%       .init.quant (def centers=[],clean_iters=2)
%       .noise.const/.rel (def 5e-6 / 8e-3)
%       .tv_ec=0.20, .tv_beta=1e-4, .tv_alpha=3e-4
%       .l2_alpha=6e-3
%       .maxIter=6, .maxIterInLine=20
%       .snap_final_to_levels=false, .min_patch_tri=0
%       .savePlots=true, .dpi=300
%
% Output:
%   rec struct (sigma0, sigma_rec, sigma_rec_tri, misfit, paths)

% -------- défauts sûrs --------
if nargin<4, opts = struct(); end
opts = def(opts, 'clip', [0.05, 0.50]);
opts = def(opts, 'doWarp', true);
opts = def(opts, 'savePlots', true);
opts = def(opts, 'dpi', 300);

if ~isfield(opts,'init'), opts.init = struct(); end
opts.init = def(opts.init, 'discrete', true);
opts.init = def(opts.init, 'calib_enable', true);
opts.init = def(opts.init, 'calib_clip', [0.95, 1.05]);
if ~isfield(opts.init,'quant') || ~isstruct(opts.init.quant)
    opts.init.quant = struct('centers',[],'clean_iters',2);
end

if ~isfield(opts,'noise') || ~isstruct(opts.noise)
    opts.noise = struct('const',5e-6,'rel',8e-3);
else
    opts.noise = def(opts.noise,'const',5e-6);
    opts.noise = def(opts.noise,'rel',  8e-3);
end

opts = def(opts, 'tv_ec',   0.20);
opts = def(opts, 'tv_beta', 1e-4);
opts = def(opts, 'tv_alpha',3e-4);
opts = def(opts, 'l2_alpha',6e-3);

opts = def(opts, 'maxIter',       6);
opts = def(opts, 'maxIterInLine', 20);

opts = def(opts, 'snap_final_to_levels', false);
opts = def(opts, 'min_patch_tri', 0);

% ---------- Chemins ----------
rootThis = fullfile('Outputs',patient_id,sprintf('slice_%03d',slice_this));
rootPrev = fullfile('Outputs',patient_id,sprintf('slice_%03d',slice_prev));
if ~exist(rootThis,'dir'), mkdir(rootThis); end
plotDir  = fullfile(rootThis,'plots_inverse'); if ~exist(plotDir,'dir'), mkdir(plotDir); end

meshThis = fullfile(rootThis,'mesh',sprintf('mesh_slice%03d.mat',slice_this));
meshPrev = fullfile(rootPrev,'mesh',sprintf('mesh_slice%03d.mat',slice_prev));
packThis = fullfile(rootThis,'eit_pack.mat');
packPrev = fullfile(rootPrev,'eit_pack.mat');
vizFile  = fullfile(rootThis,'viz_config.mat');

% ---------- Données forward ----------
S1 = load(packThis);                     % g,H,E,params,Imeas,sigma_tri (si sim)
g  = S1.g;  H = S1.H;  E = S1.E;
Ne = S1.params.Ne;
Umeas = S1.Imeas(:);
sigma_true_tri = []; if isfield(S1,'sigma_tri'), sigma_true_tri = S1.sigma_tri(:); end

% bornes/cmap (si dispo)
viz = struct('sigma_clim',[0.05 0.50],'cmap','turbo','displayMode','neurological');
if isfile(vizFile)
    try
        tmp = load(vizFile);
        if isfield(tmp,'viz_config'), viz = tmp.viz_config; end
    end
end

% ---------- Solveur FORWARD ----------
mesh = ForwardMesh1st(g*1e-3, H, E);     % mm -> m
fem  = EITFEM(mesh);
fem.mode = 'current';
fem.zeta = S1.params.z_contact * ones(Ne,1);
inj = EITSim.buildTrigPattern(Ne, S1.params.I_amp);
fem.Iel = inj(:);
fem.Uel = Umeas;

% Pondération WLS (bruit)
fem.SetInvGamma(opts.noise.const, opts.noise.rel);

% ---------- INIT depuis slice_prev ----------
[sigma0, initInfo] = init_from_prev_slice_discrete( ...
    meshPrev, packPrev, meshThis, packThis, struct('clip',opts.clip,'doWarp',opts.doWarp,'init',opts.init));
sigma0 = sigma0(:);
sigma0 = min(max(sigma0, opts.clip(1)), opts.clip(2));

% Calibration amplitude (optionnelle)
if istrue(opts.init.calib_enable)
    Uf = fem.SolveForwardVec(sigma0);
    a = real((Umeas(:)'*Uf(:)) / max(Uf(:)'*Uf(:),eps));
    a = min(max(a, opts.init.calib_clip(1)), opts.init.calib_clip(2));
    fprintf('  [calib init] a=%.4f\n', a);
    sigma0 = a*sigma0;
    sigma0 = min(max(sigma0, opts.clip(1)), opts.clip(2));
end

% --- diagnostique misfit/corr à l'init ---
Uf_init = fem.SolveForwardVec(sigma0);
misfit_init = norm(Uf_init - Umeas) / max(norm(Umeas), eps);
Vf0 = reshape(Uf_init, Ne, []); Vm0 = reshape(Umeas, Ne, []);
c0  = corr(Vm0(:,1)-mean(Vm0(:,1)), Vf0(:,1)-mean(Vf0(:,1)));
fprintf('[init] misfit_rel=%.4g  corr=%.3f  max|V_fwd|=%.6g V\n', ...
        misfit_init, c0, max(abs(Vf0),[],'all'));

% ---------- Priors ----------
tv_k = PriorTotalVariation(g, H, opts.tv_ec, 1, opts.tv_beta);
tv_k.alpha = opts.tv_alpha;
l2_k = PriorL2ToInit(sigma0, opts.l2_alpha);

% --- log @ init
sig = sigma0(:);
f_data0 = fem.OptimizationFunction(sig);
f_tv0   = tv_k.OptimizationFunction(sig);
f_l20   = l2_k.OptimizationFunction(sig);
fprintf('[ofun @init] data=%.3g  TV=%.3g  L2=%.3g  total=%.3g\n', ...
        f_data0, f_tv0, f_l20, f_data0+f_tv0+f_l20);

% ---------- Gauss-Newton (une passe) ----------
gn = SolverGN({fem, tv_k, l2_k});
set_ifprop(gn,'plotLinesearch',0);
set_ifprop(gn,'plotData',0);
set_ifprop(gn,'plotConvergence',0);
set_ifprop(gn,'plotIterations',1);
set_ifprop(gn,'plotUQ',0);

set_ifprop(gn,'alwaysMove',0);
set_ifprop(gn,'useLastDist',0);
set_ifprop(gn,'parabolic',0);
set_ifprop(gn,'nStepsBack',2);
set_ifprop(gn,'eStop',1e-3);
set_ifprop(gn,'maxIter',opts.maxIter);
set_ifprop(gn,'maxIterInLine',opts.maxIterInLine);

% Moniteur (facultatif)
try
    mon = IterMonitor(fem, Umeas, H, g, sigma0, sigma_true_tri);
    gn.plotter = mon;
    set_ifprop(gn,'outputFunVals',0);
catch
    mon = [];
end

% Optimisation
try
    sig = gn.Solve(sig);
catch ME
    if ~contains(ME.message,'EARLY_STOP_CLOSESTART')
        rethrow(ME);
    else
        fprintf('[early-stop] close-start condition.\n');
    end
end

% Meilleur état visité
if ~isempty(mon) && ~isempty(mon.sigma_best)
    sig = mon.sigma_best;
    fprintf('[post] best-iter misfit=%.6g retenu.\n', mon.best_misfit);
end

% ---------- Snap (optionnel) ----------
if istrue(opts.snap_final_to_levels) && exist('project_discrete','file') == 2 ...
   && isfield(initInfo,'levels') && ~isempty(initInfo.levels)
    mis_pre = norm(fem.SolveForwardVec(sig) - Umeas) / max(norm(Umeas), eps);
    sig_snap = project_discrete(H, g, sig, initInfo.levels, opts.min_patch_tri);
    mis_post = norm(fem.SolveForwardVec(sig_snap) - Umeas) / max(norm(Umeas), eps);
    if mis_post <= mis_pre*(1+5e-4)
        sig = sig_snap; fprintf('[snap] adopté (%.6g -> %.6g).\n', mis_pre, mis_post);
    else
        fprintf('[snap] rejeté (%.6g -> %.6g).\n', mis_pre, mis_post);
    end
end

% ---------- Clip non-dégradant ----------
sig_best = sig(:);
mis_best = norm(fem.SolveForwardVec(sig_best) - Umeas) / max(norm(Umeas), eps);
sig_clip = min(max(sig_best, opts.clip(1)), opts.clip(2));
mis_clip = norm(fem.SolveForwardVec(sig_clip) - Umeas) / max(norm(Umeas), eps);
if mis_clip <= mis_best + 1e-9
    sig = sig_clip; final_misfit = mis_clip;
    fprintf('[post-clip] misfit=%.6g (clip gardé)\n', final_misfit);
else
    sig = sig_best; final_misfit = mis_best;
    fprintf('[post-clip] clip dégradant (%.6g -> %.6g), état best conservé.\n', mis_best, mis_clip);
end

% ---------- Résultats / diagnostics ----------
sigma_rec      = sig(:);
sigma_rec_tri  = node2tri_avg(H, sigma_rec);
Ufwd           = fem.SolveForwardVec(sigma_rec);
misfit         = final_misfit;

outMat = fullfile(rootThis,'reconstruction_inverse.mat');
tracker = []; try, tracker = gn.tracker; end %#ok<NASGU>
save(outMat,'sigma0','sigma_rec','sigma_rec_tri','misfit','opts','tracker', ...
            'Umeas','g','H','Ne','sigma_true_tri','initInfo');

Vf  = reshape(Ufwd,Ne,[]); Vm = reshape(Umeas,Ne,[]);
vm1 = Vm(:,1)-mean(Vm(:,1)); vf1 = Vf(:,1)-mean(Vf(:,1));
fprintf('\n[Diagnostics data-fit]\n  max|V_meas| = %.6g V,  max|V_fwd(reco)| = %.6g V\n', ...
        max(abs(Vm),[],'all'), max(abs(Vf),[],'all'));
cc = 1.0; try, cc = corr(vm1,vf1); end
fprintf('  corr(V_meas, V_fwd) = %.3f\n', cc);

% ---------- Plot triptyque ----------
try
    ttl = sprintf('Reco — misfit=%.4g', misfit);
    if istrue(opts.savePlots)
        outPNG = fullfile(plotDir, 'triptyque_inverse.png');
    else
        outPNG = [];
    end
    plot_triptyque(g,H, sigma0, sigma_true_tri, sigma_rec_tri, viz, ttl, outPNG, opts.dpi);
catch ME
    warning('Plot triptyque impossible (%s).', ME.message);
end

% ---------- Sortie ----------
rec = struct('sigma0',sigma0,'sigma_rec',sigma_rec,'sigma_rec_tri',sigma_rec_tri, ...
             'data_misfit_rel',misfit,'misfit',misfit, ...
             'paths',struct('outMat',outMat,'plotDir',plotDir, ...
                            'packThis',packThis,'packPrev',packPrev));
end

% ================= helpers locaux ==================
function S = def(S,f,v)
if ~isfield(S,f) || isempty(S.(f)), S.(f) = v; end
end

function tf = istrue(v)
tf = ~isempty(v) && ((islogical(v) && v) || (isnumeric(v) && v~=0));
end

function set_ifprop(obj,prop,val)
try, if isprop(obj,prop), obj.(prop)=val; end, catch, end
end

function x = node2tri_avg(H, nod)
nod = nod(:);
x = (nod(H(:,1)) + nod(H(:,2)) + nod(H(:,3)))/3;
end

function plot_triptyque(g,H, sigma0_nod, sigma_true_tri, sigma_rec_tri, viz, titleReco, outPNG, dpi)
if nargin<9 || isempty(dpi), dpi = 300; end
if nargin<8, outPNG = []; end
if nargin<7 || isempty(titleReco), titleReco = 'Reco'; end

sigma0_tri = (sigma0_nod(H(:,1)) + sigma0_nod(H(:,2)) + sigma0_nod(H(:,3)))/3;

fh = figure('Color','w','Name','Triptyque','NumberTitle','off');
t = tiledlayout(1,3,'Padding','compact','TileSpacing','compact'); %#ok<NASGU>

% 1) Init
nexttile; hold on;
patch('Faces',H,'Vertices',g,'FaceVertexCData',sigma0_tri, ...
      'FaceColor','flat','EdgeColor','none');
triplot(H, g(:,1), g(:,2), 'Color', [0.75 0.75 0.75]);
axis equal tight; colormap(viz.cmap); caxis(viz.sigma_clim);
cb=colorbar; ylabel(cb,'S/m'); title('Init (\sigma_0)','FontWeight','bold');
apply_display_convention(gca, viz);

% 2) Cible
nexttile; hold on;
if ~isempty(sigma_true_tri)
    patch('Faces',H,'Vertices',g,'FaceVertexCData',sigma_true_tri, ...
          'FaceColor','flat','EdgeColor','none');
    triplot(H, g(:,1), g(:,2), 'Color', [0.75 0.75 0.75]);
    axis equal tight; colormap(viz.cmap); caxis(viz.sigma_clim);
    cb=colorbar; ylabel(cb,'S/m'); title('Cible (\sigma_{true})','FontWeight','bold');
    apply_display_convention(gca, viz);
else
    text(0.5,0.5,'\sigma_{true} indisponible','Units','normalized','HorizontalAlignment','center');
    axis off;
end

% 3) Reco
nexttile; hold on;
patch('Faces',H,'Vertices',g,'FaceVertexCData',sigma_rec_tri, ...
      'FaceColor','flat','EdgeColor','none');
triplot(H, g(:,1), g(:,2), 'Color', [0.75 0.75 0.75]);
axis equal tight; colormap(viz.cmap); caxis(viz.sigma_clim);
cb=colorbar; ylabel(cb,'S/m'); title(titleReco,'FontWeight','bold');
apply_display_convention(gca, viz);

if ~isempty(outPNG)
    exportgraphics(fh, outPNG, 'Resolution', dpi);
end
end

function apply_display_convention(ax, viz)
if nargin<1 || isempty(ax), ax=gca; end
mode = 'neurological';
if isstruct(viz) && isfield(viz,'displayMode'), mode = viz.displayMode; end
if strcmpi(mode,'radiological'), set(ax,'XDir','reverse'); else, set(ax,'XDir','normal'); end
end
