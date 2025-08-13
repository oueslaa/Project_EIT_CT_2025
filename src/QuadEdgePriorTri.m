classdef QuadEdgePriorTri < handle
    properties, L; lam; end
    methods
        function self = QuadEdgePriorTri(H, lam)
            self.lam = lam;
            T = size(H,1);
            E = sort([H(:,[1 2]); H(:,[2 3]); H(:,[1 3])],2);
            [E,~,ix] = unique(E,'rows');
            triOfEdge = accumarray(ix, repelem((1:T)',3,1), [], @(v){v});
            I=[]; J=[]; V=[];
            for k=1:numel(triOfEdge)
                t = triOfEdge{k};
                if numel(t)==2
                    a=t(1); b=t(2);
                    I=[I a a b b]; J=[J a b a b]; V=[V 1 -1 -1 1]; %#ok<AGROW>
                end
            end
            self.L = sparse(I,J,V,T,T); % Laplacien sur TRI
        end
        function f = OptimizationFunction(self, p)
            r = self.L*p(:); f = 0.5*self.lam*(r'*r);
        end
        function [H,G] = GetHessAndGrad(self, p)
            L = self.L; H = self.lam*(L'*L); G = H*p(:);
        end
        function Plot(~,~), end
    end
end
