classdef EITSim < handle
% EITSIM - Simulation EIT (forward uniquement) sur un maillage 2D
% -------------------------------------------------------------------------
% Objectif
%   - Simuler les mesures EIT (potentiels/voltages électrodes) pour un
%     maillage 2D, à partir d'une conductivité paramétrée **par triangle**.
%   - Convertir automatiquement sigma_TRI (P0) en sigma_NODALE (P1) via Φ.
%   - Fournir le Jacobien dU/dσ_triang analytique (via chaîne : Jnod*Φ) et
%     une version Finite Differences (FD) pour debug.
%
% Hypothèses / conventions
%   - g (sommets) et el_centers sont en **mm** (convertis en m pour le solveur).
%   - triGroup (Mt x 1) = étiquette de groupe par triangle (1=SoftTissue, 2=Heart, …).
%   - params.cond est une struct avec un champ par groupe (SoftTissue, Heart, …)
%     exprimés en S/m (unité SI).
%   - Les électrodes E proviennent d'un placement sur le bord (ex. buildElectrodesFromContour).
%
% Workflow typique
%   params = struct('Ne',16,'el_width',12,'I_amp',1e-3,'z_contact',0.01, ...
%                   'cond',struct('SoftTissue',0.2,'Heart',0.6,'Lung',0.08, ...
%                                  'Trachea',0.05,'Bone',0.02,'Other',0.2));
%   sim = EITSim(g,H,triGroup,domain,params,E,el_centers);
%   sim = sim.simulate_forward_full();
%   V = sim.Vmat;         % [Ne x K]
%   y = sim.Imeas;        % (Ne*K) x 1 (avec bruit optionnel)
%   [U,Jtri] = sim.forward_with_jacobian(sim.sigma_tri);
%
% Pièges courants
%   - Assure-toi que numel(E) == params.Ne (une cellule d'arêtes par électrode).
%   - Si cond.* n'a pas tous les champs requis, getfield(…,name) échouera.
%   - Le solveur OOEIT attend des **mètres** -> conversion mm->m faite en CTOR.
%   - Les patterns d'injection ont somme nulle et amplitude I_amp (par motif).
%
% AUCUN plot ici (batch-only). Sauvegarde via save_results().
% -------------------------------------------------------------------------

    properties
        % Géométrie / groups
        g; H; triGroup; domain; groups
        % Paramètres EIT
        params
        Ne; el_width; I_amp; z_contact; cond
        % Electrodes / solveur forward
        el_centers; E; solver
        % Données / états
        sigma_tri; sigma_used
        Imeas_clean; Imeas
        % Ajouts :
        U_cell          % cell(K,1): potentiel nodal par motif (Nnodesx1)
        Vmat            % [Ne x K] voltages électrodes par motif
        Ipat            % [Ne x K] motif d'injection
        E_mag_tri       % [Mt x 1] moyenne |E| sur motifs (placeholder)
        gradUx_tri      % idem
        gradUy_tri      % idem

        % === Nouveau : levée P0->P1 ===
        Phi             % (Nnodes x Ntri) : lissage/levée P0->P1 (moyenne pondérée aire)
    end

    %% ===================== CTOR & UTILS =====================
    methods
        function obj = EITSim(g, H, triGroup, domain, params, E, el_centers)
            % EITSIM  Constructeur : configure le solveur, les électrodes, Φ et les patterns.
            %
            % Entrées
            %   g, H, triGroup, domain : maillage et groupes
            %   params : struct avec au moins
            %       .Ne, .el_width, .I_amp, .z_contact, .cond (struct par groupe)
            %       [optionnels] .noiseRel (def 1e-3), .rngSeed (def [])
            %   E, el_centers : électrodes (cell) et centres (Ne x 2)
            %
            % Effets
            %   - crée l'objet OOEIT (EITFEM) sur un mesh en mètres
            %   - construit Φ (levée P0->P1)
            %   - prépare les patterns d'injection Ipat (Ne x K)

            if nargin < 6, E = []; end
            if nargin < 7, el_centers = []; end

            obj.g = g; obj.H = H; obj.triGroup = triGroup; obj.domain = domain;
            obj.params = params;
            obj.Ne = params.Ne; obj.el_width = params.el_width;
            obj.I_amp = params.I_amp; obj.z_contact = params.z_contact;

            assert(isfield(params,'cond') && isstruct(params.cond), ...
                'EITSim:cond -> params.cond (struct) est requis.');
            obj.cond = params.cond;

            if isfield(params,'groups')
                obj.groups = params.groups;
            else
                obj.groups = {'SoftTissue','Heart','Lung','Trachea','Bone','Other'};
            end
            if ~isfield(obj.params,'noiseRel'), obj.params.noiseRel = 1e-3; end
            if ~isfield(obj.params,'rngSeed'),  obj.params.rngSeed  = [];   end

            if ~isempty(E) && ~isempty(el_centers)
                obj.E = E; obj.el_centers = el_centers;
            else
                error('Provide E and el_centers from buildElectrodesFromContour.');
            end

            % Sanity checks basiques
            assert(size(obj.g,2)==2,'g must be Nx2');
            assert(all(obj.H(:)>=1) && all(obj.H(:)<=size(obj.g,1)),'H indices out of bounds');
            if ~isempty(obj.E), assert(numel(obj.E)==obj.Ne,'E count mismatch'); end
            if isnumeric(obj.triGroup)
                assert(numel(obj.triGroup)==size(obj.H,1),'triGroup length mismatch');
            end

            % ---------- Solveur OOEIT (inchangé) ----------
            g_m   = obj.g * 1e-3;  % mm -> m
            mesh  = ForwardMesh1st(g_m, obj.H, obj.E);
            solver = EITFEM(mesh);
            solver.zeta = obj.z_contact * ones(obj.Ne,1);
            solver.mode = 'current';
            obj.solver = solver;

            % ---------- Patterns d'injection ----------
            obj.Ipat = EITSim.buildTrigPattern(obj.Ne, obj.I_amp);
            obj.solver.Iel = obj.Ipat(:);  % vectorisé pour SolveForwardVec/SolveForward

            % ---------- Nouveau : construire Phi (levée P0->P1) ----------
            obj.Phi = EITSim.buildPhi(obj.g, obj.H);
        end

        function c = get_cond(obj,name)
            % GET_COND  Raccourci : renvoie params.cond.(name)
            assert(isfield(obj.cond,name),'Missing conductivity field');
            c = obj.cond.(name);
        end
    end

    %% ===================== FORWARD =====================
    methods
        function obj = simulate_forward_full(obj)
            % SIMULATE_FORWARD_FULL
            % 1) Construit σ_triang via triGroup + params.cond
            % 2) Lève en σ_nodal via Φ (P0 -> P1)
            % 3) Appelle SolveForwardVec pour tous les motifs d'un coup
            % 4) Ajoute du bruit relatif (noiseRel) pour produire Imeas

            % 1) sigma par triangle (P0)
            Mt = size(obj.H,1);
            st = zeros(Mt,1);
            if isnumeric(obj.triGroup)
                for k=1:numel(obj.groups)
                    name = obj.groups{k};
                    idx  = (obj.triGroup==k);
                    st(idx) = getfield(obj.cond,name); %#ok<GFLD>
                end
            else
                st(:) = obj.cond.SoftTissue;
            end
            obj.sigma_tri  = st(:);
            obj.sigma_used = obj.sigma_tri;

            % 2) sigma nodal via Phi (P0->P1)
            sigma_nod = obj.Phi * obj.sigma_tri;

            % 3) solve (TOUS les motifs d'un coup via SolveForwardVec)
            if ~isempty(obj.params.rngSeed), rng(obj.params.rngSeed); end
            K = size(obj.Ipat,2);
            obj.solver.Iel = obj.Ipat(:);           % empilement par colonnes: motif1; motif2;...
            Uv = obj.solver.SolveForwardVec(sigma_nod);
            obj.Vmat = reshape(Uv, obj.Ne, K);      % [Ne x K], colonnes = motifs
            % Ne x K
            obj.U_cell = {};                            % non utilisé ici

            % 4) Imeas (clean + bruit optionnel)
            obj.Imeas_clean = obj.Vmat(:);
            rel = max(obj.params.noiseRel, 0);
            noise = rel * max(abs(obj.Imeas_clean)) * randn(size(obj.Imeas_clean));
            obj.Imeas = obj.Imeas_clean + noise;

            % 5) Champ E (placeholders pour future extension)
            obj.E_mag_tri = []; obj.gradUx_tri = []; obj.gradUy_tri = [];
        end

        function U = SolveForwardVec(obj, sigma_tri)
            % SOLVEFORWARDVEC  Forward direct pour un σ_triang donné.
            % Entrée : sigma_tri (Mt x 1) TRIANGLES
            % Sortie : U ((Ne*K) x 1) mesures empilées
            sigma_nod = obj.Phi * sigma_tri(:);
            if isempty(obj.solver)
                error('EITSim:SolveForwardVec:NoBackend', ...
                      'Attach a FEM solver to obj.solver before calling SolveForwardVec.');
            end
            U = obj.solver.SolveForwardVec(sigma_nod);
        end

        function [U, Jtri] = forward_with_jacobian(obj, sigma_tri)
            % FORWARD_WITH_JACOBIAN  Mesures + Jacobien analytique par triangle.
            % U   : (Ne*K) x 1   — mesures empilées
            % Jtri: (Ne*K) x Mt  — dérivées dU/dsigma_tri (chaîne: Jnod * Φ)
            sigma_tri = sigma_tri(:);
            if isempty(obj.solver)
                error('EITSim:forward_with_jacobian:NoForward', ...
                     'Attach a solver exposing SolveForwardVec(σ) to compute [U,J].');
            end

            % P0 -> P1
            sigma_nod = obj.Phi * sigma_tri;

            % Résolution (met à jour l'état interne du solver si besoin)
            U = obj.solver.SolveForwardVec(sigma_nod);

            % Jacobien natif OOEIT par NOEUD (P1)
            % alreadyComputed=1 car on vient de résoudre avec sigma_nod
            Jnod = obj.solver.Jacobian(sigma_nod, 1);

            % Règle de la chaîne : dU/dθ = (dU/dσ_nod) * (dσ_nod/dθ) = Jnod * Phi
            Jtri = Jnod * obj.Phi;
        end

        function [U, Jfd] = forward_with_jacobian_fd(obj, sigma_tri)
            % FORWARD_WITH_JACOBIAN_FD  Version Finite Differences (DEBUG).
            % U   : (Ne*K) x 1
            % Jfd : (Ne*K) x Mt (Mt = #triangles)
            sigma_tri = sigma_tri(:);
            if isempty(obj.solver)
                error('EITSim:forward_with_jacobian_fd:NoForward', ...
                     'Attach a solver exposing SolveForwardVec(σ) to compute [U,J].');
            end

            % Baseline
            sigma_nod = obj.Phi * sigma_tri;
            U = obj.solver.SolveForwardVec(sigma_nod);

            m = numel(U); n = numel(sigma_tri);
            Jfd = zeros(m,n, 'like', U);
            eps_rel = 1e-3; eps_abs = 1e-4;

            for j = 1:n
                dj = max(eps_abs, eps_rel * abs(sigma_tri(j)));
                sp = sigma_tri; sp(j) = sp(j) + dj;
                Up = obj.solver.SolveForwardVec(obj.Phi * sp);
                Jfd(:,j) = (Up(:) - U(:)) / dj;
            end
        end
    end

    %% ===================== HELPERS (statics) =====================
    methods (Static)
        function injPattern = buildTrigPattern(Ne, I_amp)
            % BUILDTRIGPATTERN  Motifs d'injection trigonométriques.
            % - Chaque colonne (motif) a somme nulle (courant conservé).
            % - Normalisation max|.|=1 puis échelle I_amp.
            injPattern = zeros(Ne, Ne);
            idx = (0:Ne-1)';
            for k = 1:Ne
                phase = 2*pi*(k-1)/Ne;
                v = sin( 2*pi*idx/Ne + phase );
                v = v - mean(v);                 % somme nulle
                if max(abs(v)) > 0
                    v = v / max(abs(v));         % normalise
                end
                injPattern(:,k) = I_amp * v;
            end
        end

        function tri = node2tri_avg(H, nod)
            % NODE2TRI_AVG  Moyenne simple nodale -> élémentaire (diagnostics).
            nod = nod(:);
            tri = (nod(H(:,1)) + nod(H(:,2)) + nod(H(:,3))) / 3;
            tri = tri(:);
        end

        function nod = tri2node_avg(H, tri, N)
            % TRI2NODE_AVG  Ancienne levée moyenne (debug; non utilisée).
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

        function [gx, gy] = gradTriLinear(g, H, U)
            % GRADTRILINEAR  ∇U constant par triangle (base P1).
            U = U(:); Mt = size(H,1);
            gx = zeros(Mt,1); gy = zeros(Mt,1);
            for t=1:Mt
                v = H(t,:);
                x = g(v,1); y = g(v,2);
                A = [x(2)-x(1), x(3)-x(1); y(2)-y(1), y(3)-y(1)];
                phi = [U(v(2))-U(v(1)); U(v(3))-U(v(1))];
                gxy = A \ phi;      % ≈ ∇U
                gx(t) = gxy(1); gy(t) = gxy(2);
            end
        end

function Phi = buildPhi(g, H)
    ng = size(g,1); nH = size(H,1); gDim = size(g,2);

    % Aire/volume des éléments
    w = zeros(nH,1);
    if gDim == 2
        for e = 1:nH
            v = g(H(e,:),:);
            w(e) = 0.5*abs(det([v(2,:)-v(1,:); v(3,:)-v(1,:)]));
        end
    elseif gDim == 3
        for e = 1:nH
            v = g(H(e,:),:);
            w(e) = abs(det([v(2,:)-v(1,:); v(3,:)-v(1,:); v(4,:)-v(1,:)]))/6;
        end
    else
        error('EITSim.buildPhi: dimension %d non supportée', gDim);
    end

    % Assemblage creux
    I = zeros(3*nH,1); J = zeros(3*nH,1); V = zeros(3*nH,1);  % 2D -> 3 sommets/tri
    deg = zeros(ng,1);
    p = 0;
    k = size(H,2);
    for e = 1:nH
        nodes = H(e,:);
        for t = 1:k
            i = nodes(t);
            p = p + 1;
            I(p) = i; J(p) = e; V(p) = w(e);
            deg(i) = deg(i) + w(e);
        end
    end
    I = I(1:p); J = J(1:p); V = V(1:p);

    Phi = sparse(I, J, V, ng, nH);
    Phi = spdiags(1./max(deg, eps), 0, ng, ng) * Phi;
end

    end

    %% ===================== I/O =====================
    methods
        function save_results(obj, outdir)
            % SAVE_RESULTS  Sauvegarde un "pack" .mat minimal pour post-traitement.
            % Écrit : eit_pack.mat contenant maillage, électrodes, Vmat/Ipat,
            %         mesures avec/sans bruit, paramètres, placeholders.
            if nargin<2 || isempty(outdir), outdir = pwd; end
            if ~exist(outdir,'dir'), mkdir(outdir); end
            packFile = fullfile(outdir, 'eit_pack.mat');
            g = obj.g; H = obj.H; triGroup = obj.triGroup; %#ok<NASGU>
            E = obj.E; el_centers = obj.el_centers; %#ok<NASGU>
            sigma_tri = obj.sigma_tri; params = obj.params; domain = obj.domain; %#ok<NASGU>
            % Toujours sauver Vmat/Ipat + Imeas (provenant de Vmat(:))
            Imeas_clean = []; Imeas = [];
            if isprop(obj,'Vmat') && ~isempty(obj.Vmat)
                Vmat = obj.Vmat; %#ok<NASGU>
                Imeas_clean = obj.Vmat(:);
            else
                Vmat = []; %#ok<NASGU>
            end
            if isprop(obj,'Ipat'), Ipat = obj.Ipat; else, Ipat = []; end %#ok<NASGU>
            if isprop(obj,'Imeas') && ~isempty(obj.Imeas)
                Imeas = obj.Imeas; %#ok<NASGU>
            else
                Imeas = Imeas_clean; %#ok<NASGU>
            end
            U_cell = obj.U_cell; 
            E_mag_tri = obj.E_mag_tri; 
            gradUx_tri = obj.gradUx_tri; 
            gradUy_tri = obj.gradUy_tri;

            save(packFile, 'g','H','triGroup','E','el_centers','sigma_tri', ...
                'params','Imeas','Imeas_clean','domain', ...
                'U_cell','Vmat','Ipat','E_mag_tri','gradUx_tri','gradUy_tri');

        end
    end
end
