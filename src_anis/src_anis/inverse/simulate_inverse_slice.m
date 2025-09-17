function rec = simulate_inverse_slice(patient_this, slice_this, patient_prev, slice_prev, opts, iterDir)
% SIMULATE_INVERSE_SLICE
%   - Maillage FORWARD: pack cible (riche)
%   - Maillage INVERSE: reconstruit depuis contour cible + organes slice_prev
%   - Init σ0: posée depuis polygones (discrète, pas constante)
%   - NOSER/GN sur maillage inverse, projection sur forward pour plots/signaux

if nargin < 6, iterDir = []; end
if nargin < 5 || isempty(opts), opts = struct(); end
if ~isempty(iterDir)
    if getfield_def(opts,'cleanIterDir',false) && exist(iterDir,'dir')==7
        try, rmdir(iterDir,'s'); catch, end
    end
    if exist(iterDir,'dir')~=7, mkdir(iterDir); end
end

%% ====== Chemins ======
rootThis = fullfile('Outputs',patient_this,sprintf('slice_%03d',slice_this));
rootPrev = fullfile('Outputs',patient_prev,sprintf('slice_%03d',slice_prev));

meshFileF = fullfile(rootThis,'mesh',sprintf('mesh_slice%03d.mat',slice_this));
packFile  = fullfile(rootThis,'eit_pack.mat');
meshPrev  = fullfile(rootPrev,'mesh',sprintf('mesh_slice%03d.mat',slice_prev));

assert(isfile(packFile),  'Pack forward cible manquant: %s', packFile);
assert(isfile(meshFileF), 'Mesh forward cible manquant: %s', meshFileF);
assert(isfile(meshPrev),  'Mesh donneur manquant: %s', meshPrev);

S  = load(packFile);         % g,H,E,params,...
MF = load(meshFileF);        % g,H,triGroup,contour,...

%% ====== Defaults solveur ======
opts = set_default(opts,'param_by','node');   % 'node' ou 'tri'
opts = set_default(opts,'reg_alpha',0);
opts = set_default(opts,'iters_per_stage',0);
opts = set_default(opts,'constNoise',3e-5);
opts = set_default(opts,'relNoise',1.2e-2);
opts = set_default(opts,'snapshots_enable',true);

opts.noser_only = getfield_def(opts,'noser_only',struct()); 
opts.noser_only = set_default(opts.noser_only,'enable',false);
opts.noser_only = set_default(opts.noser_only,'iters',1);
opts.noser_only = set_default(opts.noser_only,'lambda',1e-3);
opts.noser_only = set_default(opts.noser_only,'lambda_mode','scaled-median');
opts.noser_only = set_default(opts.noser_only,'use_diag',true);
opts.noser_only = set_default(opts.noser_only,'max_rel_step',0.20);
opts.noser_only = set_default(opts.noser_only,'linesearch',false);
opts.noser_only = set_default(opts.noser_only,'alpha_l2',0);

% Clip si absent
if ~isfield(opts,'clip') || isempty(opts.clip) || ~all(isfinite(opts.clip))
    vals = struct2array(S.params.cond);
    lo = max(1e-3, min(vals)); hi = max(vals); span = hi-lo;
    opts.clip = [max(1e-3, lo-0.02*span), hi+0.02*span];
end

%% ====== Maillage FORWARD (pour plotting/projection) ======
if isfield(S,'g') && isfield(S,'H')
    gF_mm = S.g;   HF = S.H;     % mm
else
    gF_mm = MF.g;  HF = MF.H;    % mm
end
gF_m = gF_mm*1e-3;

%% ====== Maillage INVERSE (contour cible + organes prev) ======
invDir  = fullfile(rootThis,'mesh_inverse'); if ~exist(invDir,'dir'), mkdir(invDir); end
invFile = fullfile(invDir, sprintf('mesh_inverse_slice%03d_from_%s_%03d.mat', ...
                                   slice_this, patient_prev, slice_prev));

want_prev_org = getfield_def(opts,'inverse_from_prev_organs',true);
rebuild_flag  = getfield_def(opts,'rebuild_inverse_mesh',false);

need_rebuild = rebuild_flag || ~isfile(invFile);
invObj = [];
if ~need_rebuild
    tmp = load(invFile);
    has_shapes = isfield(tmp,'shapes') && isstruct(tmp.shapes);
    if want_prev_org && ~has_shapes
        need_rebuild = true;  % ancien cache sans organes
    else
        invObj = tmp;
    end
end

if want_prev_org
    if need_rebuild
        Pext = MF.contour;                        % contour ext. cible (mm)
        invObj = build_inverse_mesh_with_prev_organs( ...
                    meshPrev, Pext, slice_this, ...
                    getfield_def(opts,'inv_params',struct('targetSize',5,'Hgrad',1.4,'GeometricOrder','linear','minArea_mm2',300)), ...
                    getfield_def(opts,'tps',struct('enable',true,'lambda',1e-3)));
        invObj.prev_id   = patient_prev;
        invObj.prev_slice= slice_prev;
        save(invFile,'-struct','invObj');
    end
    gI_mm = invObj.g; HI = invObj.H;
else
    % maillage neutre sur contour seul
    iMesh = build_inverse_mesh_from_contour(MF.contour, ...
                getfield_def(opts,'inv_Hmax',5), getfield_def(opts,'inv_Hgrad',1.3));
    gI_mm = iMesh.g; HI = iMesh.H;
    invObj = struct('g',gI_mm,'H',HI,'shapes',struct(),'domain',polyshape(),'contour',iMesh.g); %#ok<NASGU>
    save(invFile,'-struct','iMesh');
end
gI_m = gI_mm*1e-3;  NnI = size(gI_m,1);

% Mapping TRI->NOEUD (utile si param_by='tri')
PhiI = buildPhi_tri2node(gI_mm, HI);

%% ====== Données/Patterns ======
% Mesures
if isfield(S,'Vmat') && ~isempty(S.Vmat)
    Umeas = S.Vmat(:);
elseif isfield(S,'Imeas') && ~isempty(S.Imeas)
    Umeas = S.Imeas(:);
elseif isfield(S,'Imeas_clean') && ~isempty(S.Imeas_clean)
    Umeas = S.Imeas_clean(:);
else
    error('Mesures introuvables dans pack.');
end
Ne = S.params.Ne;
K  = round(numel(Umeas)/Ne);
Z  = reshape(Umeas,Ne,[]); Z = Z - mean(Z,1); Umeas = Z(:);
scaleU = max(1, max(abs(Umeas))); Umeas = Umeas/scaleU;

% Patterns
if isfield(S,'Ipat') && ~isempty(S.Ipat), base = S.Ipat;
elseif exist('EITSim','class') == 8 || exist('EITSim','class') == 2
    base = EITSim.buildTrigPattern(Ne, S.params.I_amp);
else
    base = inj_pattern_trig(Ne, S.params.I_amp);
end
if size(base,2) < K, base = repmat(base,1,ceil(K/size(base,2))); end
inj = base(:,1:K);

%% ====== FEM (forward riche + inverse) ======
meshF = ForwardMesh1st(gF_m, HF, S.E);
meshF.SetInverseMesh(struct('g',gI_m,'H',HI));   % découplage inverse
fem = EITFEM(meshF); fem.mode='current';
fem.zeta = S.params.z_contact*ones(Ne,1);
fem.Iel  = inj(:); fem.Uel = Umeas(:);
fem.SetInvGamma(opts.constNoise, opts.relNoise);

%% ====== Paramétrisation θ ===============================================
param_by = lower(getfield_def(opts,'param_by','node'));
switch param_by
    case 'node'
        pushForward = @(theta) theta;
        pullBack    = @(sigma) sigma;
        J_apply     = @(Jnod) Jnod;
    case 'tri'
        pushForward = @(theta) PhiI*theta(:);                % σ = Φ θ
        pullBack    = @(sigma) node2tri_avg(HI, sigma(:));   % θ = avg(σ)
        J_apply     = @(Jnod) Jnod*PhiI;                     % Jθ
    otherwise
        error('opts.param_by must be ''node'' or ''tri''.');
end

%% ====== Init σ0 nodale (DISCRÈTE depuis polygones) ======================
order = {'SoftTissue','Heart','Lung','Trachea','Bone','Cartilage','Blood','LiverSpleenPancreas','Kidney','Other'};
default_soft = getfield_def(S.params.cond,'SoftTissue',0.35);

% vecteur numérique des niveaux, aligné sur 'order'
levels_vec = zeros(numel(order),1);
for j=1:numel(order)
    gj = order{j};
    levels_vec(j) = getfield_def(S.params.cond, gj, default_soft);
end

% labels nodaux sur le maillage INVERSE
lab_nod_I = ones(NnI,1);  % 1 = SoftTissue (background)
if want_prev_org && isfield(invObj,'shapes') && ~isempty(invObj.shapes)
    for j = 2:numel(order) % priorité à l'ordre (shapes déjà disjoints)
        gj = order{j};
        if isfield(invObj.shapes,gj) && invObj.shapes.(gj).NumRegions>0
            IN = isinterior(invObj.shapes.(gj), gI_mm(:,1), gI_mm(:,2));
            if any(IN), lab_nod_I(IN) = j; end
        end
    end
end
sigma0_inv = levels_vec(lab_nod_I);
sigma0_inv = min(max(sigma0_inv, opts.clip(1)), opts.clip(2));

% labels TRIANGULAIRES sur le maillage FORWARD (pour les plots discrets)
centF = (gF_mm(HF(:,1),:) + gF_mm(HF(:,2),:) + gF_mm(HF(:,3),:))/3;
lab_tri_F = ones(size(HF,1),1);
if want_prev_org && isfield(invObj,'shapes') && ~isempty(invObj.shapes)
    for j = 2:numel(order)
        gj = order{j};
        if isfield(invObj.shapes,gj) && invObj.shapes.(gj).NumRegions>0
            INt = isinterior(invObj.shapes.(gj), centF(:,1), centF(:,2));
            if any(INt), lab_tri_F(INt) = j; end
        end
    end
end

theta   = pullBack(sigma0_inv);
sig_inv = pushForward(theta);

%% ====== NOSER puis GN ====================================================
iter_global = 0; mu_hist = [];

% NOSER
if getfield_def(opts.noser_only,'enable',false) && getfield_def(opts.noser_only,'iters',0)>0
    nNOSER = max(1, round(opts.noser_only.iters));
    for kk=1:nNOSER
        U0  = fem.SolveForwardVec(sig_inv); %#ok<NASGU>
        Jnd = fem.Jacobian(sig_inv, 1);
        Jx  = J_apply(Jnd);
        W   = fem.InvGamma;

        lambda = getfield_def(opts.noser_only,'lambda',1e-3);
        if strcmpi(getfield_def(opts.noser_only,'lambda_mode','scaled-median'),'scaled-median')
            dS = full(diag(Jx.'*(W*Jx))); med = median(dS(dS>0));
            if isfinite(med) && med>0, lambda = lambda*med; end
        end
        ns = opts.noser_only; ns.lambda = lambda;

        [theta, sig_inv] = noser_step_theta_thetaParam(fem, theta, sig_inv, Jx, W, ns, opts.clip, pushForward);
        mu_hist(end+1,1) = lambda; %#ok<AGROW>
        iter_global = iter_global + 1;

        if getfield_def(opts,'snapshots_enable',true)
            sig_fwd = meshF.ItoF(sig_inv);
            save_iter_snapshot(iterDir, gF_mm, HF, sig_fwd, iter_global, 1);
        end
    end
end

% GN
maxIter = max(0, round(getfield_def(opts,'iters_per_stage',0)));
alpha   = max(0, getfield_def(opts,'reg_alpha',0));
for it = 1:maxIter
    U0  = fem.SolveForwardVec(sig_inv);
    r   = fem.Uel - U0;

    Jnd = fem.Jacobian(sig_inv, 1);
    Jx  = J_apply(Jnd);
    W   = fem.InvGamma;

    Sjj = Jx'*(W*Jx);
    gk  = Jx'*(W*r);
    if alpha>0
        Sjj = Sjj + 2*alpha*speye(length(theta));
        gk  = gk  - 2*alpha*theta;
    end

    dtheta = Sjj\gk;

    tau = 1.0; f0 = norm(r);
    for bt = 1:8
        theta_try = theta + tau*dtheta;
        sig_try   = pushForward(theta_try);
        sig_try   = min(max(sig_try, opts.clip(1)), opts.clip(2));
        f  = norm(fem.SolveForwardVec(sig_try) - fem.Uel);
        if f <= f0 || ~isfinite(f), break; end
        tau = 0.5*tau;
    end

    theta = theta + tau*dtheta;
    sig_inv = pushForward(theta);
    sig_inv = min(max(sig_inv, opts.clip(1)), opts.clip(2));

    iter_global = iter_global + 1;
    if getfield_def(opts,'snapshots_enable',true)
        sig_fwd = meshF.ItoF(sig_inv);
        save_iter_snapshot(iterDir, gF_mm, HF, sig_fwd, iter_global, 2);
    end
end

%% ====== Sorties / diagnostics ===========================================
sigma_rec_inv = sig_inv(:);
sigma_rec_fwd = meshF.ItoF(sigma_rec_inv);
sigma_rec_tri = node2tri_avg(HF, sigma_rec_fwd);

Ufwd = fem.SolveForwardVec(sigma_rec_inv);
Vf   = reshape(Ufwd,Ne,[]); Vf = Vf - mean(Vf,1); Ufwd = Vf(:)/scaleU;

Z  = reshape(Umeas,Ne,[]); Z = Z - mean(Z,1); UmeasN = Z(:);
misfit = norm(Ufwd - UmeasN)/max(norm(UmeasN),eps);

Vf0 = reshape(Ufwd*scaleU,Ne,[]);
Vm0 = reshape(UmeasN*scaleU,Ne,[]);
vm1 = Vm0(:,1)-mean(Vm0(:,1)); vf1 = Vf0(:,1)-mean(Vf0(:,1));
if norm(vm1)>0 && norm(vf1)>0, corr1 = corr(vm1, vf1, 'rows','complete'); else, corr1 = NaN; end
umeas_max = max(abs(Vm0),[],'all'); ufwd_max = max(abs(Vf0),[],'all');

sigma0_fwd = meshF.ItoF(sigma0_inv(:));

rec = struct();
rec.sigma0_inv     = sigma0_inv(:);
rec.sigma0         = sigma0_fwd(:);
rec.sigma_rec_inv  = sigma_rec_inv(:);
rec.sigma_rec      = sigma_rec_fwd(:);
rec.sigma_rec_tri  = sigma_rec_tri(:);
rec.misfit         = misfit;
rec.diag           = struct('umeas_max',umeas_max,'ufwd_max',ufwd_max,'corr',corr1);

% >>> initInfo bien formé pour le plot discret <<<
rec.initInfo       = struct( ...
    'init_mode', 'morpho-mesher', ...
    'levels',    levels_vec(:), ...         % VECTEUR NUMÉRIQUE
    'label_names',{order}, ...
    'lab_tri',   lab_tri_F(:), ...          % indices 1..numel(levels_vec) sur H (forward)
    'lab_nod',   [] );

rec.paths          = struct('packThis',packFile,'meshThisF',meshFileF,'meshPrev',meshPrev,'meshInv',invFile);
rec.Ne             = Ne;
rec.H              = HF;
rec.g              = gF_mm;
rec.Umeas          = Umeas(:);
rec.Ufwd           = Ufwd(:);
rec.sigma_clim     = opts.clip;
rec.scaleU         = scaleU;
rec.stage          = struct('noser_steps', numel(mu_hist), 'gn_iters', maxIter);
if ~isempty(mu_hist), rec.mu_hist = mu_hist; rec.mu = mu_hist(end); else, rec.mu = NaN; end
rec.theta_rec      = theta;
rec.ig             = gI_mm; rec.iH = HI;
end

%% ===================== Helpers =====================
function S = set_default(S,f,v), if ~isfield(S,f) || isempty(S.(f)), S.(f)=v; end, end
function v = getfield_def(S,f,def), if isstruct(S)&&isfield(S,f)&&~isempty(S.(f)), v=S.(f); else, v=def; end, end

function tri = node2tri_avg(H, nod)
nod = nod(:);
tri = (nod(H(:,1))+nod(H(:,2))+nod(H(:,3)))/3;
end

function Phi = buildPhi_tri2node(g, H)
ng = size(g,1); Mt = size(H,1);
w = zeros(Mt,1);
for e=1:Mt
    v=g(H(e,:),:);
    w(e)=0.5*abs(det([v(2,:)-v(1,:); v(3,:)-v(1,:)]));
end
I=zeros(3*Mt,1); J=I; V=I; deg=zeros(ng,1); p=0;
for e=1:Mt
    nodes=H(e,:);
    for t=1:3
        i=nodes(t); p=p+1;
        I(p)=i; J(p)=e; V(p)=w(e);
        deg(i)=deg(i)+w(e);
    end
end
I=I(1:p); J=J(1:p); V=V(1:p);
Phi=sparse(I,J,V,ng,Mt);
deg=max(deg,eps);
Phi=spdiags(1./deg,0,ng,ng)*Phi;
end

function save_iter_snapshot(iterDir,g_mm,H,sig_nodal_fwd,iter_global,stage_k)
if isempty(iterDir) || exist(iterDir,'dir')~=7, return; end
tri = (sig_nodal_fwd(H(:,1))+sig_nodal_fwd(H(:,2))+sig_nodal_fwd(H(:,3)))/3;
S=struct('sigma_iter_nodal',sig_nodal_fwd(:),'sigma_iter_tri',tri(:), ...
         'iter_global',iter_global,'stage_k',stage_k);
save(fullfile(iterDir,sprintf('iter_%03d.mat',iter_global)),'-struct','S');
end

function [theta_new, sig_new] = noser_step_theta_thetaParam(fem, theta, sig, Jx, W, ns, clip, pushForward)
r = fem.Uel - fem.SolveForwardVec(sig);
S = Jx.'*(W*Jx);
if getfield_def(ns,'use_diag',true), R=spdiags(diag(S),0,size(S,1),size(S,2));
else, R=speye(size(S,1)); end
alpha = getfield_def(ns,'alpha_l2',0);
H = S + getfield_def(ns,'lambda',1e-3)*R + 2*alpha*speye(size(S,1));
g = Jx.'*(W*r) - 2*alpha*theta;

dtheta = H\g;
mrs = getfield_def(ns,'max_rel_step',0.10);
t_norm=norm(theta); d_norm=norm(dtheta);
if isfinite(mrs) && d_norm>mrs*max(t_norm,1e-12) && d_norm>0
    dtheta = dtheta*((mrs*max(t_norm,1e-12))/d_norm);
end

if getfield_def(ns,'linesearch',false)
    tau=1.0; f0=norm(r);
    for it=1:6
        theta_try = theta + tau*dtheta;
        sig_try   = pushForward(theta_try);
        sig_try   = min(max(sig_try, clip(1)), clip(2));
        f  = norm(fem.SolveForwardVec(sig_try) - fem.Uel);
        if f<=f0 || ~isfinite(f), break; end
        tau=0.5*tau;
    end
    theta_new = theta + tau*dtheta;
else
    theta_new = theta + dtheta;
end
sig_new = pushForward(theta_new);
sig_new = min(max(sig_new, clip(1)), clip(2));
end

function injPattern = inj_pattern_trig(Ne, I_amp)
injPattern = zeros(Ne, Ne);
idx = (0:Ne-1)';
for k=1:Ne
    phase = 2*pi*(k-1)/Ne;
    v = sin(2*pi*idx/Ne + phase);
    v = v - mean(v);
    if max(abs(v))>0, v = v/max(abs(v)); end
    v = v - mean(v);
    injPattern(:,k) = I_amp*v;
end
end

function iMesh = build_inverse_mesh_from_contour(P_ext, Hmax, Hgrad)
if nargin<2 || isempty(Hmax),  Hmax  = 5; end
if nargin<3 || isempty(Hgrad), Hgrad = 1.3; end
P = clean_open_poly(P_ext);
k=size(P,1);
gd=zeros(2+2*k,1); gd(1)=2; gd(2)=k;
gd(3:2+k)=P(:,1); gd(3+k:2+2*k)=P(:,2);
ns=char('P1')'; sf='P1';
model=createpde(); [dl,~]=decsg(gd,sf,ns); geometryFromEdges(model,dl);
msh=generateMesh(model,'Hmax',Hmax,'Hgrad',Hgrad,'GeometricOrder','linear','Jiggle',{'on','minimum'});
iMesh.g=msh.Nodes.'; iMesh.H=msh.Elements.';
end 

function P = clean_open_poly(Pin)
P = Pin;
tol = 1e-6;
d=[inf; hypot(diff(P(:,1)), diff(P(:,2)))]; P(d<tol,:)=[];
P=unique(P,'rows','stable'); if ~isequal(P(1,:),P(end,:)), P(end+1,:)=P(1,:); end; P(end,:)=[];
end
