classdef IterMonitor < handle
    properties
        fem; Umeas; H; sigma0; sigmaTrueTri = []
        iter = 0
        last_misfit = NaN
        best_misfit = Inf
        sigma_best  = []
        last_sigma  = []      % pour la dérive incrémentale

        drift_cap   = Inf     % cumul vers σ0 (mettre Inf si on ne veut pas l'utiliser)
        step_cap    = 0.12    % dérive relative par pas: ||σ_k-σ_{k-1}||/||σ_{k-1}||
        print_corr  = true
    end
    methods
        function self = IterMonitor(fem,Umeas,H,~,sigma0,sigmaTrueTri,varargin)
            self.fem=fem; self.Umeas=Umeas(:); self.H=H; self.sigma0=sigma0(:);
            if ~isempty(sigmaTrueTri), self.sigmaTrueTri=sigmaTrueTri(:); end
            for k=1:2:numel(varargin)
                p=varargin{k}; v=varargin{k+1}; if isprop(self,p), self.(p)=v; end
            end
        end
        function plot(self,sigma)
            self.iter = self.iter + 1; sigma = sigma(:);

            Uf = self.fem.SolveForwardVec(sigma);
            misfit = norm(Uf - self.Umeas)/max(norm(self.Umeas),eps);

            rel_to_init = norm(sigma - self.sigma0)/max(norm(self.sigma0),eps);
            if isempty(self.last_sigma), step_rel = 0;
            else, step_rel = norm(sigma - self.last_sigma)/max(norm(self.last_sigma),eps);
            end

            s_corr = NaN;
            if self.print_corr
                try
                    Ne = self.fem.fmesh.nEl;
                    Vf = reshape(Uf,Ne,[]); Vm = reshape(self.Umeas,Ne,[]);
                    s_corr = corr(Vm(:,1)-mean(Vm(:,1)), Vf(:,1)-mean(Vf(:,1)));
                catch, s_corr = NaN; end
            end

            relErr_true = NaN;
            if ~isempty(self.sigmaTrueTri)
                tri = (sigma(self.H(:,1))+sigma(self.H(:,2))+sigma(self.H(:,3)))/3;
                relErr_true = norm(tri - self.sigmaTrueTri)/max(norm(self.sigmaTrueTri),eps);
            end

            gain_per_drift = NaN;
            if ~isnan(self.last_misfit)
                dm = self.last_misfit - misfit;
                gain_per_drift = dm / max(step_rel,1e-12);
            end
            self.last_misfit = misfit;

            if misfit < self.best_misfit, self.best_misfit=misfit; self.sigma_best=sigma; end

            if isnan(s_corr)
                fprintf('[iter %d] misfit=%.6g  step=%.3g  ||Δσ||/||σ0||=%.3g  gain/step=%.3g', ...
                    self.iter, misfit, step_rel, rel_to_init, gain_per_drift);
            else
                fprintf('[iter %d] misfit=%.6g  corr=%.3f  step=%.3g  ||Δσ||/||σ0||=%.3g  gain/step=%.3g', ...
                    self.iter, misfit, s_corr, step_rel, rel_to_init, gain_per_drift);
            end
            if ~isnan(relErr_true), fprintf('  relErr_true=%.3g', relErr_true); end
            fprintf('\n');

            % arrêt seulement si le pas est trop gros ET n'apporte rien
            if (step_rel > self.step_cap) && (~isnan(gain_per_drift) && gain_per_drift < 1e-5)
                error('IterMonitor:EarlyStop','EARLY_STOP_CLOSESTART');
            end
            % (drift_cap cumulatif optionnel)
            if ~isinf(self.drift_cap) && (rel_to_init > self.drift_cap)
                error('IterMonitor:EarlyStop','EARLY_STOP_CLOSESTART');
            end

            self.last_sigma = sigma;
        end
    end
end
