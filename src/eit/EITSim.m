classdef EITSim < handle
% EITSIM - Simulation EIT (forward uniquement) sur un maillage 2D
% - Paramétrisation interne: sigma par TRIANGLE (un triangle = un organe).
% - Le solveur FEM (EITFEM) attend sigma NODALE -> conversion TRI->NOD côté EITSim.
% - AUCUN plot ici (batch-only).
%
% Dépendances externes :
%   - ForwardMesh1st.m
%   - EITFEM.m  (expose SolveForwardVec(sigma_nod))

    properties
        % Géométrie / groups
        g; H; triGroup; domain; groups
        % Paramètres EIT
        params
        Ne; el_width; I_amp; z_contact; cond
        % Electrodes / solveur forward
        el_centers; E; solver
        % Données / états
        sigma_tri;       % (nH x 1) sigma par élément (triangle)
        sigma_used
        Imeas_clean; Imeas
        sigma_rec
        sigma_init_used   % [Nnodes x 1] S/m (optionnel)
    end

    %% ===================== CTOR & UTILS =====================
    methods
        function obj = EITSim(g, H, triGroup, domain, params, E, el_centers)
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

            % Sanity checks
            assert(size(obj.g,2)==2,'g must be Nx2');
            assert(all(obj.H(:)>=1) && all(obj.H(:)<=size(obj.g,1)),'H indices out of bounds');
            if ~isempty(obj.E), assert(numel(obj.E)==obj.Ne,'E count mismatch'); end

            if isnumeric(obj.triGroup)
                assert(numel(obj.triGroup)==size(obj.H,1),'triGroup length mismatch');
            end
        end

        function c = get_cond(obj,name)
            assert(isfield(obj.cond,name),'Missing conductivity field');
            c = obj.cond.(name);
        end
    end

    %% ===================== FORWARD (backend requis) =====================
    methods
        function obj = simulate_forward(obj)
            % 1) sigma_tri par triangle à partir de triGroup/cond
            Mtri = size(obj.H,1);
            st = zeros(Mtri,1);

            if isnumeric(obj.triGroup)
                K = numel(obj.groups);
                for k=1:K
                    name = obj.groups{k};
                    idx  = (obj.triGroup==k);
                    if isfield(obj.cond,name), st(idx) = obj.cond.(name);
                    else, st(idx) = obj.cond.SoftTissue; end
                end
            elseif isstruct(obj.triGroup)
                flds = fieldnames(obj.triGroup);
                st(:) = obj.cond.SoftTissue;
                for i=1:numel(flds)
                    name = flds{i};
                    tris = obj.triGroup.(name);
                    tris = tris(tris>=1 & tris<=Mtri);
                    if isfield(obj.cond,name), st(tris) = obj.cond.(name); end
                end
            else
                st(:) = obj.cond.SoftTissue;
            end
            obj.sigma_tri  = st(:);      % (Mt x 1) par élément (organe)
            obj.sigma_used = obj.sigma_tri;

            % 2) backend FEM (init si besoin, sinon réutilise)
            if isempty(obj.solver)
                g_m   = obj.g * 1e-3;                 % mm -> m pour le FEM
                mesh  = ForwardMesh1st(g_m, obj.H, obj.E);
                solver = EITFEM(mesh);
                solver.zeta = obj.z_contact * ones(obj.Ne,1);
                solver.mode = 'current';
                injPattern  = EITSim.buildTrigPattern(obj.Ne, obj.I_amp);
                solver.Iel  = injPattern(:);
                obj.solver = solver;
            else
                if ~isprop(obj.solver,'Iel') || isempty(obj.solver.Iel)
                    injPattern = EITSim.buildTrigPattern(obj.Ne, obj.I_amp);
                    obj.solver.Iel = injPattern(:);
                end
                if isprop(obj.solver,'zeta') && (isempty(obj.solver.zeta) || numel(obj.solver.zeta)~=obj.Ne)
                    obj.solver.zeta = obj.z_contact * ones(obj.Ne,1);
                end
            end

            % 3) Conversion ELEMENT -> NŒUD pour le solveur FEM (EITFEM)
            sigma_nod = EITSim.tri2node_avg(obj.H, obj.sigma_tri, size(obj.g,1));

            % 4) Mesures
            obj.Imeas_clean = obj.solver.SolveForwardVec(sigma_nod);  % EITFEM (inchangé)
            if ~isempty(obj.params.rngSeed), rng(obj.params.rngSeed); end
            rel = max(obj.params.noiseRel, 0);
            noise = rel * max(abs(obj.Imeas_clean)) * randn(size(obj.Imeas_clean));
            obj.Imeas = obj.Imeas_clean + noise;
        end

        function U = SolveForwardVec(obj, sigma_tri)
            % Wrapper pour forward direct à sigma par TRIANGLE (Mt x 1).
            % Convertit en sigma NODALE et délègue à EITFEM.SolveForwardVec.
            sigma_nod = EITSim.tri2node_avg(obj.H, sigma_tri(:), size(obj.g,1));
            if isempty(obj.solver)
                error('EITSim:SolveForwardVec:NoBackend', ...
                      'Attach a FEM solver to obj.solver before calling SolveForwardVec.');
            end
            U = obj.solver.SolveForwardVec(sigma_nod);  % EITFEM (inchangé)
        end

        function [U, J] = forward_with_jacobian(obj, sigma_tri)
            % Fallback Jacobien par différences finies en param TRIANGLE.
            % (utile si le backend ne fournit pas J ; ne modifie pas EITFEM)
            sigma_tri = sigma_tri(:);
            if isempty(obj.solver)
                error('EITSim:forward_with_jacobian:NoForward', ...
                     'Attach a solver exposing SolveForwardVec(σ) to compute [U,J].');
            end
            sigma_nod = EITSim.tri2node_avg(obj.H, sigma_tri, size(obj.g,1));
            U = obj.solver.SolveForwardVec(sigma_nod);

            m = numel(U); n = numel(sigma_tri);
            J = zeros(m,n, 'like', U);
            eps_rel = 1e-3; eps_abs = 1e-4;
            for j = 1:n
                dj = max(eps_abs, eps_rel * abs(sigma_tri(j)));
                sp = sigma_tri; sp(j) = sp(j) + dj;
                Up = obj.solver.SolveForwardVec( ...
                        EITSim.tri2node_avg(obj.H, sp, size(obj.g,1)));
                J(:,j) = (Up - U) / dj;
            end
        end
    end

    %% ===================== HELPERS (statics) =====================
    methods (Static)
        function injPattern = buildTrigPattern(Ne, I_amp)
            injPattern = zeros(Ne, Ne);
            for k = 1:Ne
                phase = 2*pi*(k-1)/Ne;
                for n = 1:Ne
                    injPattern(n,k) = I_amp * sin( 2*pi*(n-1)/Ne + phase );
                end
            end
        end

        function tri = node2tri_avg(H, nod)
            nod = nod(:);
            tri = (nod(H(:,1)) + nod(H(:,2)) + nod(H(:,3))) / 3;
            tri = tri(:);
        end

        function nod = tri2node_avg(H, tri, N)
            % Moyenne des valeurs élémentaires sur les nœuds adjacents.
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

        function Ctri = sigma_tri_from_groups(H, triGroup, groups, cond)
            M = size(H,1);
            Ctri = zeros(M,1);
            for k = 1:numel(groups)
                name = groups{k};
                idx  = (triGroup == k);
                if isfield(cond,name)
                    Ctri(idx) = cond.(name);
                else
                    Ctri(idx) = cond.SoftTissue;
                end
            end
        end
    end

    %% ===================== I/O (pack-only, pas de plots) =====================
    methods
        function save_results(obj, outdir)
            % Sauvegarde un pack .mat standard **sans aucun PNG**.
            if nargin<2 || isempty(outdir), outdir = pwd; end
            if ~exist(outdir,'dir'), mkdir(outdir); end

            packFile = fullfile(outdir, 'eit_pack.mat');
            g = obj.g; H = obj.H; triGroup = obj.triGroup;
            E = obj.E; el_centers = obj.el_centers;
            sigma_tri = obj.sigma_tri; params = obj.params; domain = obj.domain; %#ok<NASGU>
            Imeas_clean = []; Imeas = [];
            if isprop(obj,'Imeas_clean'), Imeas_clean = obj.Imeas_clean; end
            if isprop(obj,'Imeas'),       Imeas       = obj.Imeas;       end

            save(packFile, 'g','H','triGroup','E','el_centers','sigma_tri', ...
                            'params','Imeas','Imeas_clean','domain');
            fprintf('[EITSim] Pack sauvegardé (mat-only): %s\n', packFile);
        end
    end
end
