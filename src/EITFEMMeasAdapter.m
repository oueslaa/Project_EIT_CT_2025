classdef EITFEMMeasAdapter < handle
    % Emballe un objet EITFEM pour lui fournir l'écart Data–Model
    properties
        base    % l'objet EITFEM d'OOEIT
        meas    % Imeas (dimension compatible avec Uel)
        W       % matrice de poids (Ln' InvGamma Ln) + diag loading
        name = 'EITFEM(data misfit)';
    end
    methods
        function self = EITFEMMeasAdapter(base, meas)
            self.base = base;
            self.meas = meas(:);
            % récupérer Ln et InvGamma du EITFEM
            Ln = base.Ln;            % [m x m]
            G  = base.InvGamma;      % [m x m]
            self.W = Ln' * G * Ln;
            % --- diagonal loading pour stabilité ---
            self.W = self.W + 1e-8 * speye(size(self.W,1));
        end
        function val = OptimizationFunction(self, sigma)
            U = self.base.SolveForwardVec(sigma);  % Uel
            r = U - self.meas;                     % résidu
            val = 0.5 * (r'*(self.W*r));
        end
        function [Hess, grad] = GetHessAndGrad(self, sigma)
            % Jacobien J (potentiels électrode vs sigma)
            J = self.base.Jacobian(sigma);
            % Hess ≈ J' W J
            Hess = J' * (self.W * J);
            % grad = J' W (U - meas)
            U = self.base.SolveForwardVec(sigma);
            grad = J' * (self.W * (U - self.meas));
        end
        function Plot(self, sigma)
            % Reutilise le plot de EITFEM (données vs modèle)
            self.base.Plot(sigma);
        end
    end
end
