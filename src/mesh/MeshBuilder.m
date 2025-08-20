classdef MeshBuilder
    properties
        g; H; triGroup; shapes; domain; groups; group_names; contour;
    end

    methods(Static)
        function obj = from_polys(shapes, domain, ext_smooth, z_slice, params)
            % Params par défaut
            if ~isfield(params, 'targetSize'),   params.targetSize   = 5;   end
            if ~isfield(params, 'minArea_mm2'),  params.minArea_mm2  = 300; end
            if ~isfield(params, 'Hgrad'),        params.Hgrad        = 1.4; end
            if ~isfield(params, 'GeometricOrder'), params.GeometricOrder = 'linear'; end

            types = {'SoftTissue','Heart','Lung','Trachea','Bone','Other'};
            names = {'Soft Tissue','Heart','Lung','Trachea','Bone','Other'};
            tol = 1e-6;

            % -------- 1) Contour extérieur nettoyé --------
            P_ext = ext_smooth;
            dP    = [inf; hypot(diff(P_ext(:,1)), diff(P_ext(:,2)))]; 
            P_ext(dP<tol,:) = [];
            P_ext = unique(P_ext, 'rows', 'stable');
            if ~isequal(P_ext(1,:), P_ext(end,:)), P_ext(end+1,:) = P_ext(1,:); end
            P_ext(end,:) = [];
            if size(P_ext,1) < 3, error('Contour extérieur insuffisant pour maillage.'); end

            cleanContours = {P_ext}; 
            groupList = {'SoftTissue'};

            % -------- 2) Organes par régions --------
            orgFields = fieldnames(shapes);
            for k = 1:numel(orgFields)
                shp = shapes.(orgFields{k}); 
                if shp.NumRegions == 0, continue; end
                regs = regions(shp); grp = orgFields{k};
                for r = 1:numel(regs)
                    P = regs(r).Vertices;
                    dP = [inf; hypot(diff(P(:,1)), diff(P(:,2)))]; P(dP<tol,:) = [];
                    P = unique(P, 'rows', 'stable');
                    if ~isequal(P(1,:), P(end,:)), P(end+1,:) = P(1,:); end
                    P(end,:) = [];
                    if size(P,1) < 3 || abs(polyarea(P(:,1),P(:,2))) < params.minArea_mm2, continue; end
                    cleanContours{end+1} = P; %#ok<AGROW>
                    groupList{end+1}     = grp; %#ok<AGROW>
                end
            end

            % -------- 3) Géométrie PDE --------
            nS = numel(cleanContours); 
            vk = cellfun(@(C) size(C,1), cleanContours); 
            maxk = max(vk);
            gd = zeros(2 + 2*maxk, nS);
            for i = 1:nS
                P = cleanContours{i}; k = size(P,1);
                gd(1,i) = 2; gd(2,i) = k; 
                gd(3:2+k,i)       = P(:,1); 
                gd(3+k:2+2*k,i)   = P(:,2);
            end
            sf = strjoin(arrayfun(@(i) sprintf('P%d', i), 1:nS, 'UniformOutput', false), '+');
            ns = char(arrayfun(@(i) sprintf('P%d', i), 1:nS, 'UniformOutput', false))';

            model = createpde();
            [dl, ~] = decsg(gd, sf, ns); 
            geometryFromEdges(model, dl);

            % -------- 4) Maillage --------
            msh = generateMesh(model, ...
                'Hmax', params.targetSize, ...
                'Hgrad', params.Hgrad, ...
                'GeometricOrder', params.GeometricOrder, ...
                'Jiggle', {'on','minimum'});
            g = msh.Nodes'; 
            H = msh.Elements'; 
            contour = P_ext;

            % -------- 5) triGroup vectorisé par barycentres --------
            M = size(H,1);
            centroids = ( g(H(:,1),:)+g(H(:,2),:)+g(H(:,3),:) )/3;
            triGroup = ones(M,1); % 1 = SoftTissue par défaut

            % Affectation vectorisée : pour chaque organe interne (2..end)
            for k = 2:numel(types)
                if isfield(shapes, types{k})
                    IN = isinterior(shapes.(types{k}), centroids(:,1), centroids(:,2)); % logique vectorisée
                    triGroup(IN) = k;
                end
            end

            % Sécurité : si un centroid n'est pas dans "domain", repasse en SoftTissue
            try
                OUT = ~isinterior(domain, centroids(:,1), centroids(:,2));
                triGroup(OUT) = 1;
            catch
                % si 'domain' n'est pas polyshape (rare), on ignore
            end

            % -------- 6) Objet --------
            obj = MeshBuilder; 
            obj.g = g; obj.H = H; obj.triGroup = triGroup; 
            obj.shapes = shapes; obj.domain = domain; 
            obj.groups = types; obj.group_names = names; 
            obj.contour = contour;
        end
    end 

    methods
        function saveMesh(obj, outdir, z_slice)
            outdir = char(outdir);
            if ~exist(outdir, 'dir'), mkdir(outdir); end
            fname = fullfile(outdir, sprintf('mesh_slice%03d.mat', z_slice));
            g = obj.g; H = obj.H; triGroup = obj.triGroup; 
            shapes = obj.shapes; domain = obj.domain; contour = obj.contour;
            save(fname, 'g', 'H', 'triGroup', 'shapes', 'domain', 'contour');
            fprintf('[MeshBuilder] Maillage sauvegardé : %s\n', fname);
        end

        function show(obj)
            figure('Color','w'); hold on;
            cmap = lines(numel(obj.groups));
            for k = 1:numel(obj.groups)
                idx = (obj.triGroup == k); 
                if ~any(idx), continue; end
                tris = obj.H(idx,:);
                patch('Faces', tris, 'Vertices', obj.g, ...
                    'FaceColor', cmap(k,:), 'FaceAlpha', 0.85, 'EdgeColor', 'none');
            end
            P = obj.contour; 
            if ~isequal(P(1,:), P(end,:)), P=[P;P(1,:)]; end
            plot(P(:,1), P(:,2), 'r-', 'LineWidth', 1.5);
            axis equal tight; 
            title('Mesh et groupes'); 
            legend(obj.group_names, 'Location','BestOutside');
        end
    end
end
