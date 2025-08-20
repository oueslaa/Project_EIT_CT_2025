classdef PriorL2ToInit < handle
% Prior quadratique: (lambda/2) * || W .* (sigma - sigma0) ||^2
% Compatible SolverGN (OOEIT): expose OptimizationFunction() et GetHessAndGrad()

    properties
        sigma0   % vecteur nodal de référence (ng x 1)
        W        % poids (ng x 1) ou scalaire
        lambda   % force globale (scalaire)
        ng
    end

    methods
        function obj = PriorL2ToInit(sigma0, lambda, W)
            if nargin < 3 || isempty(W), W = 1; end
            obj.sigma0 = sigma0(:);
            obj.lambda = lambda;
            obj.W = W;
            obj.ng = numel(obj.sigma0);
        end

        function val = OptimizationFunction(self, sigma)
            d = sigma(:) - self.sigma0;
            if isscalar(self.W)
                val = 0.5 * self.lambda * (self.W^2) * (d.'*d);
            else
                w = self.W(:);
                val = 0.5 * self.lambda * sum( (w .* d).^2 );
            end
        end

        function [Hess, grad] = GetHessAndGrad(self, sigma)
            d = sigma(:) - self.sigma0;
            if isscalar(self.W)
                grad = self.lambda * (self.W^2) * d;
                % Hessien diagonal = lambda * W^2 * I
                Hess = self.lambda * (self.W^2) * speye(self.ng);
            else
                w = self.W(:);
                grad = self.lambda * (w.^2) .* d;
                Hess = spdiags(self.lambda * (w.^2), 0, self.ng, self.ng);
            end
        end
    end
end
