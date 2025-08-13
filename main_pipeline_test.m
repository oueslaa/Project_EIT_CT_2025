%% main_test.m — Reconstruction slice 301 en réutilisant les résultats enregistrés
clear; close all; clc;
rmpath(genpath('src'));                                  % (optionnel) repartir propre
addpath(genpath('/Users/anis/Documents/StageInria/Code/OOEIT-main'), '-begin');
addpath('src','-begin'); addpath(genpath('src'),'-begin');
rehash; clear classes;

% ===== Réglages =====
set(groot,'defaultFigureVisible','on');                  % 'off' si run silencieux
displayMode = 'neurological';
PLOT = struct('dpi',300,'pdf',true,'svg',false,'fig',false,'mode',displayMode);

patient_id = 's0011';
tgt_slice  = 301;                                        % slice à reconstruire
src_slice  = 302;                                        % slice utilisée pour sigma_init

rootTgt = fullfile('Outputs', patient_id, sprintf('slice_%03d', tgt_slice));
rootSrc = fullfile('Outputs', patient_id, sprintf('slice_%03d', src_slice));

mesh_mat_tgt = fullfile(rootTgt, 'mesh', sprintf('mesh_slice%03d.mat', tgt_slice));
res_mat_tgt  = fullfile(rootTgt, 'eit_results.mat');     % mesures + E + params de 301
res_mat_src  = fullfile(rootSrc, 'eit_results.mat');     % carte forward / résultats de 302

assert(exist(mesh_mat_tgt,'file')==2, 'Mesh cible introuvable: %s', mesh_mat_tgt);
assert(exist(res_mat_src,'file')==2,  'Résultats source introuvables: %s', res_mat_src);

% ===== 1) Charger mesh cible =====
M = load(mesh_mat_tgt);                                  % g, H, triGroup, domain
g = M.g; H = M.H; triGroup = M.triGroup; domain = M.domain;
groups = {'SoftTissue','Heart','Lung','Trachea','Bone','Other'};

% ===== 2) Params / mesures / électrodes =====
E = []; el_centers = []; params = struct(); Imeas = [];
if exist(res_mat_tgt,'file')
    T = load(res_mat_tgt);
    if isfield(T,'E'),           E = T.E; end
    if isfield(T,'el_centers'),  el_centers = T.el_centers; end
    if isfield(T,'Imeas'),       Imeas = T.Imeas; end
    params = struct();
    params.Ne        = getfield_or(T,'Ne',16);
    params.el_width  = getfield_or(T,'el_width',33);
    params.I_amp     = getfield_or(T,'I_amp',0.35e-3);
    params.z_contact = getfield_or(T,'z_contact',0.05);
    params.groups    = groups;
    params.noiseRel  = getfield_or(T,'params',struct(),'noiseRel',1e-3);
    params.rngSeed   = getfield_or(T,'params',struct(),'rngSeed',[]);
    if isfield(T,'cond')
        params.cond = T.cond;
    else
        params.cond = struct('SoftTissue',0.3,'Heart',0.5,'Lung',0.15,'Trachea',0.15,'Bone',0.05,'Other',0.4);
    end
else
    warning('eit_results (cible) manquant: %s. Paramètres par défaut.', res_mat_tgt);
    params = struct('Ne',16,'el_width',33,'I_amp',0.35e-3,'z_contact',0.05, ...
                    'groups',{groups}, 'cond',struct('SoftTissue',0.3,'Heart',0.5,'Lung',0.15,'Trachea',0.15,'Bone',0.05,'Other',0.4), ...
                    'noiseRel',1e-3,'rngSeed',[]);
end

% Si le bruit réel est inconnu, on fait un peu plus confiance au prior
params.noiseRel = max(params.noiseRel, 5e-3);

if isempty(E) || isempty(el_centers)
    [E, el_centers] = buildElectrodesFromContour(g, [], params.Ne, params.el_width, H);
end

% ===== 3) σ_init = projection carte FORWARD (σ_tri) de 302 =====
lb = 0.05; ub = 0.50;                                     % bornes S/m
sigma0 = sigma_init_from_forward_slice(res_mat_src, g, domain, params.cond.SoftTissue, [lb ub]);

% lissage léger + recentrage doux (ta "meilleure version")
sigma0 = smooth_nodal_on_mesh(H, sigma0, 2, 0.20);
sigma0 = 0.90*sigma0 + 0.10*params.cond.SoftTissue;

% >>> Nettoyage des bords d'organes : on remplace la bordure par la valeur interne
[sigma0, repSnap] = clean_sigma_borders(H, sigma0, triGroup, ...
    'ring', 1, ...          % halo de 1-anneau autour des frontières
    'method', 'median', ... % valeur interne = médiane
    'blend', 0.0, ...       % 0 = snap dur; 0.3 = mélange doux
    'bounds', [lb ub]);
fprintf('sigma0 (from forward slice %d) stats: min=%.4f  max=%.4f  std=%.4f  | border nodes snapped=%d\n', ...
        src_slice, min(sigma0), max(sigma0), std(sigma0), repSnap.nBorder);

% ===== 4) Solveur FEM (identique au modèle des mesures) =====
eitSim = EITSim(g, H, triGroup, domain, params, E, el_centers);

g_m   = g * 1e-3;
mesh  = ForwardMesh1st(g_m, H, E);
solver = EITFEM(mesh);
solver.zeta = params.z_contact * ones(params.Ne,1);
solver.mode = 'current';
inj = EITSim.buildTrigPattern(params.Ne, params.I_amp);
solver.Iel = inj(:);
eitSim.solver = solver;

if ~isempty(Imeas)
    eitSim.Imeas = Imeas(:);
else
    warning('Imeas manquant: forward homogène (SoftTissue).');
    eitSim.sigma_tri   = params.cond.SoftTissue * ones(size(H,1),1);
    eitSim.Imeas_clean = solver.SolveForwardVec(eitSim.sigma_tri);
    eitSim.Imeas       = eitSim.Imeas_clean;
end

% ===== 5) Inverse (figures GN activées) =====
% tik_weight un peu renforcé pour coller au prior nettoyé
reg = struct('lb',lb,'ub',ub,'tik_weight',2.0);
eitSim = eitSim.simulate_inverse(sigma0, true, reg);

% Clip final
eitSim.sigma_rec = min(max(eitSim.sigma_rec, lb), ub);

% ===== 6) Figures "overview" =====
expDir = fullfile(rootTgt, sprintf('from_slice_%03d', src_slice), 'overview');
if ~exist(expDir,'dir'), mkdir(expDir); end
eitSim.plot_all(expDir, displayMode, PLOT);

% ===== 7) Exports dédiés =====
expDir = fullfile(rootTgt, sprintf('from_slice_%03d', src_slice));
if ~exist(expDir,'dir'), mkdir(expDir); end

clim_init = prctile(sigma0,[1 99]); if diff(clim_init)==0, clim_init=[min(sigma0) max(sigma0)]; end
rec_vals  = eitSim.sigma_rec(:);
clim_rec  = prctile(rec_vals,[1 99]); if diff(clim_rec)==0, clim_rec=[min(rec_vals) max(rec_vals)]; end
delta     = eitSim.sigma_rec - sigma0;
L         = max(abs(prctile(delta,[1 99]))); if ~isfinite(L) || L==0, L=0.05; end
clim_delta = [-L L];

% σ_init
fig = figure('Color','w'); ax=axes('Parent',fig); hold(ax,'on');
draw_field(ax,H,g,sigma0); axis(ax,'equal'); axis(ax,'tight'); colormap(ax,jet); colorbar(ax);
caxis(ax,clim_init); title(ax,'\sigma_{init} (inverse)','Interpreter','tex');
try, apply_display_convention(ax,displayMode); end
save_plot_simple(fig, fullfile(expDir,'sigma_init'), PLOT);

% σ_rec
fig = figure('Color','w'); ax=axes('Parent',fig); hold(ax,'on');
draw_field(ax,H,g,eitSim.sigma_rec); axis(ax,'equal'); axis(ax,'tight'); colormap(ax,jet); colorbar(ax);
caxis(ax,clim_rec); title(ax,'\sigma reconstruite','Interpreter','tex');
try, apply_display_convention(ax,displayMode); end
save_plot_simple(fig, fullfile(expDir,'sigma_rec'), PLOT);

% Δσ
fig = figure('Color','w'); ax=axes('Parent',fig); hold(ax,'on');
draw_field(ax,H,g,delta); axis(ax,'equal'); axis(ax,'tight'); colormap(ax,parula); colorbar(ax);
caxis(ax,clim_delta); title(ax,'\Delta\sigma = \sigma_{rec} - \sigma_{init}','Interpreter','tex');
try, apply_display_convention(ax,displayMode); end
save_plot_simple(fig, fullfile(expDir,'delta_sigma'), PLOT);

% Sauvegarde
S = struct();
S.g=eitSim.g; S.H=eitSim.H; S.triGroup=eitSim.triGroup; S.domain=eitSim.domain; S.params=eitSim.params;
S.Ne=eitSim.Ne; S.el_width=eitSim.el_width; S.I_amp=eitSim.I_amp; S.z_contact=eitSim.z_contact;
S.cond=eitSim.cond; S.groups=eitSim.groups; S.el_centers=eitSim.el_centers; S.E=eitSim.E;
S.sigma_tri=eitSim.sigma_tri; S.sigma_used=eitSim.sigma_used;
S.Imeas_clean=eitSim.Imeas_clean; S.Imeas=eitSim.Imeas;
S.sigma_rec=eitSim.sigma_rec; S.sigma_init_used = sigma0;
save(fullfile(expDir,'eit_results_from_src.mat'),'-struct','S','-v7.3');

disp('main_test terminé : reconstruction 301 avec sigma_init = carte forward de 302');

%% ===================== Helpers locaux =====================

function val = getfield_or(S, f, defStruct, subf, def)
if nargin==3
    if isfield(S,f), val = S.(f); else, val = defStruct; end
else
    if isfield(S,f) && isfield(S.(f),subf), val = S.(f).(subf); else, val = def; end
end
end

function save_plot_simple(fig, out_noext, opts)
if ~exist(fileparts(out_noext),'dir'), mkdir(fileparts(out_noext)); end
exportgraphics(fig, [out_noext '.png'], 'Resolution', opts.dpi, 'BackgroundColor','white');
if isfield(opts,'pdf') && opts.pdf
    exportgraphics(fig, [out_noext '.pdf'], 'ContentType','image', 'BackgroundColor','white');
end
if isfield(opts,'fig') && opts.fig, savefig(fig, [out_noext '.fig']); end
end

function draw_field(ax, H, g, C)
N = size(g,1); M = size(H,1); C = C(:);
if numel(C)==N
    patch('Faces',H,'Vertices',g,'FaceVertexCData',C,'FaceColor','interp', ...
          'EdgeColor',[0.75 0.75 0.75],'EdgeAlpha',0.3,'Parent',ax);
elseif numel(C)==M
    patch('Faces',H,'Vertices',g,'FaceVertexCData',C,'FaceColor','flat', ...
          'EdgeColor',[0.75 0.75 0.75],'EdgeAlpha',0.3,'Parent',ax);
else
    Cnod = tri2node_avg(H,C,N);
    patch('Faces',H,'Vertices',g,'FaceVertexCData',Cnod,'FaceColor','interp', ...
          'EdgeColor',[0.75 0.75 0.75],'EdgeAlpha',0.3,'Parent',ax);
end
end

function nodal = tri2node_avg(H, triVals, N)
triVals = triVals(:);
if numel(triVals) ~= size(H,1)
    if mod(numel(triVals), size(H,1)) == 0
        triVals = mean(reshape(triVals, [], size(H,1))', 2);
    else
        triVals = triVals(1:size(H,1));
    end
end
nodal = zeros(N,1); cnt = zeros(N,1);
for t = 1:size(H,1)
    v = H(t,:);
    nodal(v) = nodal(v) + triVals(t); cnt(v) = cnt(v) + 1;
end
idx = cnt>0; nodal(idx) = nodal(idx)./cnt(idx);
end

function sigma0 = sigma_init_from_forward_slice(res_mat_src, g_tgt, domain_tgt, soft_val, bounds)
if nargin<5 || isempty(bounds), bounds=[0,inf]; end
lb=bounds(1); ub=bounds(2);
S = load(res_mat_src);
if     isfield(S,'sigma_tri'),  sig_src = S.sigma_tri;
elseif isfield(S,'sigma_used'), sig_src = S.sigma_used;
elseif isfield(S,'sigma_rec'),  sig_src = S.sigma_rec;
else, warning('Aucune carte trouvée. σ_init homogène.'); sigma0 = soft_val*ones(size(g_tgt,1),1); return;
end
if ~(isfield(S,'g') && isfield(S,'H'))
    warning('Résultats source sans g/H. σ_init homogène.'); sigma0 = soft_val*ones(size(g_tgt,1),1); return;
end
g_src=S.g; H_src=S.H;
if numel(sig_src)==size(H_src,1)
    sig_src = tri2node_avg(H_src, sig_src, size(g_src,1));
elseif numel(sig_src)~=size(g_src,1)
    warning('Taille inattendue. σ_init homogène.');
    sigma0=soft_val*ones(size(g_tgt,1),1); return;
end
try
    F = scatteredInterpolant(g_src(:,1), g_src(:,2), sig_src, 'linear','nearest');
catch
    F = scatteredInterpolant(g_src(:,1), g_src(:,2), sig_src, 'nearest','nearest');
end
sigma0 = F(g_tgt(:,1), g_tgt(:,2));
bad = ~isfinite(sigma0); if any(bad), sigma0(bad)=soft_val; end
try
    in = isinterior(domain_tgt, g_tgt(:,1), g_tgt(:,2));
    sigma0(~in)=soft_val;
end
sigma0 = min(max(sigma0,lb),ub);
end

function nodeOrg = node_majority_group(H, triGroup, N)
% organe majoritaire par nœud à partir des triangles adjacents
adj = cell(N,1);
for t=1:size(H,1)
    v = H(t,:); adj{v(1)}(end+1)=t; adj{v(2)}(end+1)=t; adj{v(3)}(end+1)=t;
end
nodeOrg = ones(N,1);
for i=1:N
    tg = triGroup(adj{i});
    if ~isempty(tg), nodeOrg(i) = mode(tg); end
end
end

function [sigmaOut, rep] = clean_sigma_borders(H, sigmaIn, triGroup, varargin)
% Remplace la bordure des organes par la valeur interne (médiane/ moyenne).
%   - 'ring'   : # d'anneaux autour de la frontière (1 recommandé)
%   - 'method' : 'median' (defaut) ou 'mean'
%   - 'blend'  : 0 (snap dur) .. 1 (ne rien changer)
%   - 'bounds' : [lb ub] clamp final
p = inputParser; p.KeepUnmatched=true;
addParameter(p,'ring',1);
addParameter(p,'method','median');
addParameter(p,'blend',0.0);
addParameter(p,'bounds',[0 inf]);
parse(p,varargin{:});
ring = max(0, round(p.Results.ring));
useMedian = strcmpi(p.Results.method,'median');
blend = max(0,min(1,p.Results.blend));
lb = p.Results.bounds(1); ub = p.Results.bounds(2);

N = max(H(:));
sigma = sigmaIn(:);

% --- voisinage noeuds + groupes majoritaires
E = sort([H(:,[1 2]); H(:,[2 3]); H(:,[1 3])], 2);
E = unique(E,'rows');
neigh = accumarray([E(:) [E(:,2);E(:,1)]], 1, [N N], @sum, 0)>0; % bool adjacency
nodeOrg = node_majority_group(H, triGroup, N);

% --- noeuds de frontière : arêtes où l'organe change
isBorder = false(N,1);
chg = nodeOrg(E(:,1)) ~= nodeOrg(E(:,2));
if any(chg)
    bnodes = unique(E(chg,:));
    isBorder(bnodes) = true;
end

% --- dilatation en anneaux si demandé
if ring > 0
    cur = isBorder;
    for r=1:ring
        touch = (neigh*double(cur))>0;
        isBorder = isBorder | touch;
        cur = touch;
    end
end

% --- valeur interne par organe (sur noeuds non-frontière)
sigmaOut = sigma;
K = max(nodeOrg);
for k=1:K
    maskOrg    = (nodeOrg==k);
    interior   = maskOrg & ~isBorder;
    if ~any(maskOrg), continue; end
    if ~any(interior), interior = maskOrg; end    % fallback si organe trop petit
    val = useMedian*median(sigma(interior)) + (~useMedian)*mean(sigma(interior));
    idx = maskOrg & isBorder;
    if any(idx)
        sigmaOut(idx) = (1-blend)*val + blend*sigma(idx);
    end
end

% --- clamp
sigmaOut = min(max(sigmaOut,lb),ub);
rep = struct('nBorder', nnz(isBorder));
end
