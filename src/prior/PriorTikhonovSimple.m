classdef PriorTikhonovSimple < handle
    % J(σ) = 1/2 * λ * (σ-σref)' * (L+εI) * (σ-σref)
    properties
        L
        lambda
        sigmaRef
    end
    methods
        function obj = PriorTikhonovSimple(H, N, lambda, sigmaRef, epsd)
            if nargin < 5, epsd = 1e-6; end
            obj.lambda  = lambda;
            obj.sigmaRef = sigmaRef(:);
            % Laplacien de graphe (adjacence des noeuds)
            edges = [H(:,[1 2]); H(:,[2 3]); H(:,[3 1])];
            A = sparse(edges(:,1), edges(:,2), 1, N, N);
            A = max(A, A'); % symétriser
            d = sum(A,2);
            L = spdiags(d,0,N,N) - A;
            obj.L = L + epsd*speye(N); % SPD
        end
        function val = OptimizationFunction(obj, sigma)
            r = sigma(:) - obj.sigmaRef;
            val = 0.5 * obj.lambda * (r'*(obj.L*r));
        end
        function [Hess, grad] = GetHessAndGrad(obj, sigma)
            r = sigma(:) - obj.sigmaRef;
            Hess = obj.lambda * obj.L;
            grad = obj.lambda * (obj.L*r);
        end
        function Plot(~, ~) % rien à tracer
        end
    end
end
