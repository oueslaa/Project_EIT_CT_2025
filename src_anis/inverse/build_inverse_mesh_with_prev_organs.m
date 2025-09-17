function out = build_inverse_mesh_with_prev_organs(meshPrevFile, ext_target_mm, slice_this, params, tps)
% BUILD_INVERSE_MESH_WITH_PREV_ORGANS
%  - Warpe les organes (polyshape) de la slice précédente dans le contour
%    extérieur de la slice cible, les rend disjoints, les clippe au contour
%    cible, puis remaillage via MeshBuilder.from_polys.
%
% out fields: g,H,shapes,domain,contour

if nargin<5 || isempty(tps), tps = struct('enable',true,'lambda',1e-3); end
if nargin<4 || isempty(params)
    params = struct('targetSize',5,'Hgrad',1.4,'GeometricOrder','linear','minArea_mm2',300);
end

Mprev = load(meshPrevFile);                           % doit contenir shapes, contour (mm)
assert(isfield(Mprev,'shapes'),'meshPrevFile ne contient pas "shapes".');
assert(~isempty(ext_target_mm) && size(ext_target_mm,2)==2, 'ext_target_mm invalide.');

% --- Nettoyage contours (ouverts) ----------------------------------------
Pext  = clean_open_poly(ext_target_mm);              % contour cible (301)
Pprev = []; 
if isfield(Mprev,'contour') && ~isempty(Mprev.contour)
    Pprev = clean_open_poly(Mprev.contour);          % contour donneur (310)
end

% --- Warp des organes slice_prev -> contour cible ------------------------
shPrev = rmholes_all(Mprev.shapes);
if ~isempty(Pprev)
    if getfield_def(tps,'enable',true)
        M = min(max(200, size(Pprev,1)), 800);
        Psrc = resample_closed(Pprev, M);
        Ptgt = resample_closed(Pext,  M);
        model = fit_tps(Psrc, Ptgt, getfield_def(tps,'lambda',1e-3));
        shWarp = warp_shapes_tps(shPrev, model);
    else
        T = fit_similarity(Pprev, Pext);
        shWarp = warp_shapes_affine(shPrev, T);
    end
else
    warning('Contour prev absent, pas de warp.');
    shWarp = shPrev;
end

% --- Clip dans le contour cible + disjonction par priorité ----------------
dom = polyshape([Pext; Pext(1,:)], 'Simplify', true);

shI = intersect_each(shWarp, dom);
order = {'Heart','Lung','Trachea','Bone','Cartilage','Blood','LiverSpleenPancreas','Kidney','Other'};
shI = make_disjoint_in_order(shI, order);
shI = rmholes_all(shI);
shI = prune_small_islands(shI, getfield_def(params,'minArea_mm2',300));

% --- Domaine soft-tissue = ext - union(organes) --------------------------
U = polyshape();
flds = fieldnames(shI);
for i=1:numel(flds), U = union(U, shI.(flds{i})); end
domain = subtract(dom, U); domain = simplify(domain,'KeepCollinear',true);

% --- Remaillage inverse ---------------------------------------------------
mb = MeshBuilder.from_polys(shI, domain, [Pext; Pext(1,:)], slice_this, params);

out = struct('g',mb.g,'H',mb.H,'shapes',shI,'domain',domain,'contour',mb.contour);
end

%% ================= helpers géométrie =================
function P = clean_open_poly(Pin)
P = Pin;
tol = 1e-6;
d=[inf; hypot(diff(P(:,1)), diff(P(:,2)))]; P(d<tol,:)=[];
P=unique(P,'rows','stable');
if ~isequal(P(1,:),P(end,:)), P(end+1,:)=P(1,:); end
P(end,:)=[];
end

function Q = resample_closed(P, M)
if ~isequal(P(1,:),P(end,:)), P=[P;P(1,:)]; end
s = [0; cumsum(vecnorm(diff(P),2,2))];
L = s(end); su = linspace(0,L,M).';
Q = interp1(s, P, su, 'linear');
end

function model = fit_tps(P, Q, lambda)
N=size(P,1); K = tpsK(P,P); K(1:N+1:end)=K(1:N+1:end)+lambda;
P1=[P,ones(N,1)];
A=[K,P1; P1', zeros(3,3)];
Bx=[Q(:,1);0;0;0]; By=[Q(:,2);0;0;0];
sx=A\Bx; sy=A\By;
model=struct('P',P,'wX',sx(1:N),'ax',sx(N+1:N+3),'wY',sy(1:N),'ay',sy(N+1:N+3));
end

function K=tpsK(X,Y)
NX=size(X,1); NY=size(Y,1); K=zeros(NX,NY);
for i=1:NX
    d=hypot(X(i,1)-Y(:,1), X(i,2)-Y(:,2)); d(d==0)=1;
    v=(d.^2).*log(d); v(~isfinite(v))=0; K(i,:)=v;
end
if NX==NY, ii=1:NX; K(sub2ind([NX,NY],ii,ii))=0; end
end

function Y=tps_apply(model,X)
M=size(X,1); K=tpsK(X,model.P); X1=[X,ones(M,1)];
Y=[K*model.wX + X1*model.ax,  K*model.wY + X1*model.ay];
end

function shW = warp_shapes_tps(sh, model)
flds=fieldnames(sh); shW=struct();
for i=1:numel(flds)
    P=sh.(flds{i});
    if P.NumRegions==0, shW.(flds{i})=P; continue; end
    regs=regions(P); U=polyshape();
    for r=1:numel(regs)
        V=regs(r).Vertices; if isempty(V), continue; end
        if ~isequal(V(1,:),V(end,:)), V=[V;V(1,:)]; end
        V=unique(V,'rows','stable'); if size(V,1)<3, continue; end
        W=tps_apply(model,V); W=unique(W,'rows','stable');
        Q=polyshape(W,'Simplify',true);
        U=union(U,Q);
    end
    shW.(flds{i})=simplify(U,'KeepCollinear',true);
end
end

function T = fit_similarity(P,Q)
Pc=mean(P,1); Qc=mean(Q,1); X=P-Pc; Y=Q-Qc;
[U,~,V]=svd(X.'*Y); R=U*V.'; if det(R)<0, U(:,end)=-U(:,end); R=U*V.'; end
s = trace((R.'*(X.'*Y))) / max(trace(X.'*X), eps);
if ~isfinite(s)||s<=0, s=1; end
t = Qc.' - s*R*Pc.';
T = eye(3); T(1:2,1:2)=s*R; T(1:2,3)=t;
end

function shW = warp_shapes_affine(sh,T)
flds=fieldnames(sh); shW=struct();
for i=1:numel(flds)
    P=sh.(flds{i}); if P.NumRegions==0, shW.(flds{i})=P; continue; end
    regs=regions(P); U=polyshape();
    for r=1:numel(regs)
        V=regs(r).Vertices; if isempty(V), continue; end
        if ~isequal(V(1,:),V(end,:)), V=[V;V(1,:)]; end
        V=unique(V,'rows','stable'); if size(V,1)<3, continue; end
        Vh=[V,ones(size(V,1),1)]*T.'; Vh=unique(Vh,'rows','stable');
        Q=polyshape(Vh(:,1:2),'Simplify',true);
        U=union(U,Q);
    end
    shW.(flds{i})=simplify(U,'KeepCollinear',true);
end
end

function shI = intersect_each(sh, dom)
flds=fieldnames(sh); shI=struct();
for i=1:numel(flds)
    P=sh.(flds{i});
    if P.NumRegions==0, shI.(flds{i})=polyshape(); continue; end
    Q=intersect(P, dom);
    shI.(flds{i}) = simplify(Q,'KeepCollinear',true);
end
end

function shD = make_disjoint_in_order(sh, order)
shD=sh; U=polyshape();
for j=1:numel(order)
    gj=order{j};
    if ~isfield(shD,gj), shD.(gj)=polyshape(); continue; end
    P=shD.(gj); if P.NumRegions==0, continue; end
    P2=subtract(P,U); P2=simplify(P2,'KeepCollinear',true);
    shD.(gj)=P2;
    U=union(U,P2);
end
end

function sh = rmholes_all(sh)
flds=fieldnames(sh);
for i=1:numel(flds)
    sh.(flds{i}) = rmholes_single(sh.(flds{i}));
end
end

function P2 = rmholes_single(P1)
if P1.NumRegions==0, P2=P1; return; end
regs=regions(P1); P2=polyshape();
for r=1:numel(regs)
    V=regs(r).Vertices;
    if isempty(V), continue; end
    isn=any(isnan(V),2);
    if any(isn), V=V(1:find(isn,1,'first')-1,:); end
    if size(V,1)<3, continue; end
    if ~isequal(V(1,:),V(end,:)), V=[V;V(1,:)]; end
    Q=polyshape(V,'Simplify',true);
    P2=union(P2,Q);
end
P2=simplify(P2,'KeepCollinear',true);
end

function sh = prune_small_islands(sh, Amin)
if Amin<=0, return; end
flds=fieldnames(sh);
for i=1:numel(flds)
    P=sh.(flds{i});
    if P.NumRegions==0, continue; end
    regs=regions(P);
    Q=polyshape();
    for r=1:numel(regs)
        A=abs(polyarea(regs(r).Vertices(:,1), regs(r).Vertices(:,2)));
        if A>=Amin
            Q=union(Q, rmholes_single(reconstruct_polyshape(regs(r))));
        end
    end
    sh.(flds{i}) = simplify(Q,'KeepCollinear',true);
end
end

function P = reconstruct_polyshape(R)
V=R.Vertices;
if isempty(V), P=polyshape(); return; end
isn=any(isnan(V),2);
if any(isn), V=V(1:find(isn,1,'first')-1,:); end
if size(V,1)<3, P=polyshape(); return; end
if ~isequal(V(1,:),V(end,:)), V=[V;V(1,:)]; end
P = polyshape(V,'Simplify',true);
end

function v = getfield_def(S,f,def), if isstruct(S)&&isfield(S,f)&&~isempty(S.(f)), v=S.(f); else, v=def; end, end
