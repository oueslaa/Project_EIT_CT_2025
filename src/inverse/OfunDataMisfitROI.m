classdef OfunDataMisfitROI < handle
% OfunDataMisfitROI: fonctionnelle de données pour SolverGN
% - Paramétrisation par région (une valeur de sigma par groupe d'organes)
% - Ne dépend PAS d'EITFEM.OptimizationFunction/GetHessAndGrad
%   => on utilise seulement solver.SolveForwardVec(sigma_tri)
%
% f(p) = 1/2 || L ( U(p) - y ) ||^2 + 1/2 * lambda || W^(1/2) (p - p0) ||^2
%
% p(k) = conductivité associée au groupe k
% U(p) = tensions aux électrodes simulées par le forward
%
% Requis par SolverGN:
%   - OptimizationFunction(p)
%   - GetHessAndGrad(p)  (GN approx via différences finies sur p (K DOF))
% Optionnel:
%   - Plot(~) (non utilisée ici)

    properties
        solver              % handle EITFEM
        H                   % triangles [Mt x 3]
        groupIdx            % [Mt x 1] entier dans 1..K (groupe de chaque triangle)
        y                   % mesure alignée [mU x 1]
        L                   % pondération (mU x mU) (sparse)
        p0                  % prior mean [K x 1]
        W                   % prior metric (K x K) (sparse diagonal)
        lambda              % poids Tikhonov
        epsFD               % pas pour diff. finies (relatif)
        lastU               % cache dernier U(p) (optionnel)
        lastp               % cache dernier p (optionnel)
    end

    methods
        function self = OfunDataMisfitROI(solver, H, groupIdx, y, L, reg)
            self.solver   = solver;
            self.H        = H;
            self.groupIdx = groupIdx(:);
            self.y        = y(:);
            if nargin < 5 || isempty(L), L = speye(numel(self.y)); end
            self.L        = L;

            K = max(self.groupIdx);
            % régularisation
            if nargin < 6 || isempty(reg), reg = struct(); end
            if ~isfield(reg,'tik_weight'), reg.tik_weight = 0; end
            if ~isfield(reg,'p0')         , reg.p0         = zeros(K,1); end
            if ~isfield(reg,'W')          , reg.W          = speye(K); end

            self.lambda = reg.tik_weight;
            self.p0     = reg.p0(:);
            if isscalar(reg.W), self.W = speye(K)*reg.W; else, self.W = reg.W; end
            self.epsFD  = 1e-2; % 1% relatif par défaut
            self.lastU  = [];
            self.lastp  = [];
        end

        function f = OptimizationFunction(self, p)
            U = self.forward_from_p(p);
            r = U - self.y;
            f_data = 0.5 * sum((self.L * r).^2);
            dp = p(:) - self.p0;
            f_reg  = 0.5 * self.lambda * (dp' * (self.W * dp));
            f = f_data + f_reg;
        end

        function [Hess, grad] = GetHessAndGrad(self, p)
            % Gauss-Newton approx avec J (mU x K) par différences finies (centrées)
            U0 = self.forward_from_p(p);
            r  = U0 - self.y;
            mU = numel(U0);
            K  = max(self.groupIdx);
            J  = zeros(mU, K);
            p  = p(:);

            % Pas relatif par composante
            epsk = self.epsFD * max(1, abs(p));

            for k = 1:K
                pkp = p; pkm = p;
                pkp(k) = p(k) + epsk(k);
                pkm(k) = p(k) - epsk(k);
                Up = self.forward_from_p(pkp);
                Um = self.forward_from_p(pkm);
                J(:,k) = (Up - Um) / (2*epsk(k));
            end

            Lr = self.L' * (self.L * r);
            grad = J.' * Lr + self.lambda * (self.W * (p - self.p0));

            % Hessien GN: J' L' L J + lambda W
            JL = self.L * J;
            Hess = JL.' * JL + self.lambda * self.W;
        end

        % ---- Helpers ----
        function U = forward_from_p(self, p)
            % map p -> sigma_tri -> U
            sig_tri = p(self.groupIdx); % [Mt x 1]
            % EITFEM attend un sigma par triangle (comme SolveForwardVec dans ton code)
            U = self.solver.SolveForwardVec(sig_tri(:));
            self.lastU = U; self.lastp = p(:);
        end
    end
end
