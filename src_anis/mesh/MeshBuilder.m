classdef MeshBuilder
% MESHBUILDER
% -------------------------------------------------------------------------
% Objet utilitaire pour construire un maillage 2D (PDE Toolbox) d'une coupe
% CT à partir :
%   - du contour extérieur lissé (ext_smooth),
%   - des polyshapes d'organes par groupes (shapes),
%   - du domaine "soft tissue" (domain = extérieur - organes),
% puis :
%   - générer la géométrie PDE (decsg),
%   - trianguler (generateMesh),
%   - labelliser chaque triangle suivant le groupe anatomique (triGroup),
%   - sauvegarder/afficher le résultat.
%
% Hypothèses & unités
% -------------------
% - Toutes les coordonnées passées (ext_smooth, shapes, domain) sont en mm,
%   dans la même convention d'affichage (radiologique/neurologique).
% - Les polyshapes de `shapes` sont déjà "cohérents" (peu de recouvrement)
%   et sont supposés contenus dans le contour extérieur.
%
% Workflow typique
% ----------------
% [domain, shapes, ext, info] = extract_slice_polys(patient, z, opts);
% mb = MeshBuilder.from_polys(shapes, domain, ext, z, struct('targetSize',5));
% mb.saveMesh(fullfile(rootOut,'mesh'), z);
% mb.show();
%
% Pièges courants
% ---------------
% - Contour extérieur trop court (<3 points) -> erreur "insuffisant".
% - Groupes vides : c'est OK, ils sont ignorés.
% - Triangles dont le centroïde sort de `domain` : reclassés en SoftTissue.
% - Paramètre `targetSize` trop grand -> maillage trop grossier.
%
% Propriétés exposées
% -------------------
% g         : [Ng x 2] sommets (x,y) en mm
% H         : [Mt x 3] connectivité triangles (indices dans g)
% triGroup  : [Mt x 1] étiquette de groupe (1=SoftTissue, 2=Heart, …)
% shapes    : struct de polyshape (copie d'entrée)
% domain    : polyshape du "soft tissue" (copie d'entrée)
% groups    : cellstr codes internes des groupes (ordre = codes triGroup)
% group_names: labels lisibles (même ordre que groups)
% contour   : [N x 2] contour extérieur utilisé (mm, ouvert – sans doublon final)
%
% ©2025

    properties
        g; H; triGroup;     % maillage: sommets, éléments, étiquettes triangle
        shapes; domain;     % géométrie "logique" (polyshapes)
        groups; group_names;% noms des groupes (codes & labels)
        contour;            % contour extérieur utilisé pour bâtir la géométrie
    end

    methods (Static)
        function obj = from_polys(shapes, domain, ext_smooth, z_slice, params)
            % FROM_POLYS  Construit un MeshBuilder à partir de polyshapes/contours.
            %
            % Entrées
            % -------
            % shapes     : struct de polyshape par groupe (champs attendus ci-dessous)
            % domain     : polyshape du "soft tissue" (extérieur - organes)
            % ext_smooth : [N x 2] contour extérieur lissé (mm), fermé (1er==dernier)
            % z_slice    : indice de la coupe (stocké à titre informatif si besoin)
            % params     : struct optionnelle :
            %    .targetSize      (def 5 mm)   -> Hmax de generateMesh (taille max tri)
            %    .minArea_mm2     (def 300)    -> filtre régions trop petites
            %    .Hgrad           (def 1.4)    -> gradient de taille des éléments
            %    .GeometricOrder  (def 'linear') -> 'linear' ou 'quadratic'
            %
            % Groupes reconnus
            % ----------------
            % types = {'SoftTissue','Heart','Lung','Trachea','Bone','Cartilage',...
            %          'Blood','LiverSpleenPancreas','Kidney','Other'};
            %
            % Sortie
            % ------
            % obj : instance MeshBuilder initialisée (g,H,triGroup,meta).
            %
            % Notes
            % -----
            % - La géométrie PDE (decsg) est construite à partir d'une liste
            %   de polygones (contour ext. + régions d'organes).
            % - L'affectation des triangles aux groupes est vectorisée via
            %   les centroïdes et `isinterior(polyshape, x, y)`.

            % Params par défaut
            if ~isfield(params, 'targetSize'),   params.targetSize   = 5;   end
            if ~isfield(params, 'minArea_mm2'),  params.minArea_mm2  = 300; end
            if ~isfield(params, 'Hgrad'),        params.Hgrad        = 1.4; end
            if ~isfield(params, 'GeometricOrder'), params.GeometricOrder = 'linear'; end

            types = {'SoftTissue','Heart','Lung','Trachea','Bone','Cartilage','Blood','LiverSpleenPancreas','Kidney','Other'};
            names = {'Soft Tissue','Heart','Lung','Trachea','Bone','Cartilage','Blood','Liver/Spleen/Pancreas','Kidney','Other'};

            tol = 1e-6;

            % -------- 1) Contour extérieur nettoyé --------------------------------
            % - enlève segments quasi nuls & doublons exacts,
            % - impose un contour "ouvert" pour decsg (dernier != premier).
            P_ext = ext_smooth;
            dP    = [inf; hypot(diff(P_ext(:,1)), diff(P_ext(:,2)))]; 
            P_ext(dP<tol,:) = [];
            P_ext = unique(P_ext, 'rows', 'stable');
            if ~isequal(P_ext(1,:), P_ext(end,:)), P_ext(end+1,:) = P_ext(1,:); end
            P_ext(end,:) = [];
            if size(P_ext,1) < 3, error('Contour extérieur insuffisant pour maillage.'); end

            cleanContours = {P_ext}; 
            groupList = {'SoftTissue'};

            % -------- 2) Organes par régions --------------------------------------
            % Convertit chaque région de chaque polyshape en polygone (ouvert)
            % et filtre les toutes petites surfaces (< minArea_mm2).
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

            % -------- 3) Géométrie PDE (decsg / geometryFromEdges) ----------------
            % Assemble un "Geometry Description" (gd) avec des polygones ouverts.
            nS = numel(cleanContours); 
            vk = cellfun(@(C) size(C,1), cleanContours); 
            maxk = max(vk);
            gd = zeros(2 + 2*maxk, nS);
            for i = 1:nS
                P = cleanContours{i}; k = size(P,1);
                gd(1,i) = 2; gd(2,i) = k;      % 2 = polygone
                gd(3:2+k,i)       = P(:,1);    % x
                gd(3+k:2+2*k,i)   = P(:,2);    % y
            end
            % Chaque polygone est appelé P1, P2, ... et on les "additionne" :
            sf = strjoin(arrayfun(@(i) sprintf('P%d', i), 1:nS, 'UniformOutput', false), '+');
            ns = char(arrayfun(@(i) sprintf('P%d', i), 1:nS, 'UniformOutput', false))';

            model = createpde();
            [dl, ~] = decsg(gd, sf, ns); 
            geometryFromEdges(model, dl);

            % -------- 4) Maillage --------------------------------------------------
            % Hmax contrôle la taille max des triangles (mm). Hgrad lisse la variation.
            msh = generateMesh(model, ...
                'Hmax', params.targetSize, ...
                'Hgrad', params.Hgrad, ...
                'GeometricOrder', params.GeometricOrder, ...
                'Jiggle', {'on','minimum'});
            g = msh.Nodes'; 
            H = msh.Elements'; 
            contour = P_ext;

            % -------- 5) triGroup (étiquetage par centroïdes) ---------------------
            % Vectorisé : on calcule les centroïdes de tous les triangles
            % puis on teste l'appartenance à chaque polyshape de groupe.
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

            % Sécurité : si un centroïde sort de `domain`, on re-classe en SoftTissue
            try
                OUT = ~isinterior(domain, centroids(:,1), centroids(:,2));
                triGroup(OUT) = 1;
            catch
                % si 'domain' n'est pas polyshape (cas improbable), on ignore
            end

            % -------- 6) Construction de l'objet ----------------------------------
            obj = MeshBuilder; 
            obj.g = g; obj.H = H; obj.triGroup = triGroup; 
            obj.shapes = shapes; obj.domain = domain; 
            obj.groups = types; obj.group_names = names; 
            obj.contour = contour;
        end
    end 

    methods
        function saveMesh(obj, outdir, z_slice)
            % SAVEMESH  Sauvegarde le maillage et ses méta-données sur disque.
            %
            % Entrées
            % -------
            % outdir  : dossier de sortie (créé si nécessaire)
            % z_slice : indice de la coupe (utilisé dans le nom de fichier)
            %
            % Format écrit
            % -----------
            % mesh_slice%03d.mat contenant :
            %   - g, H, triGroup, shapes, domain, contour
            %
            outdir = char(outdir);
            if ~exist(outdir, 'dir'), mkdir(outdir); end
            fname = fullfile(outdir, sprintf('mesh_slice%03d.mat', z_slice));
            g = obj.g; H = obj.H; triGroup = obj.triGroup; 
            shapes = obj.shapes; domain = obj.domain; contour = obj.contour;
            p = fileparts(fname);
            if ~isfolder(p), [ok,msg] = mkdir(p); assert(ok, "mkdir('%s') a échoué: %s", p, msg); end

            save(fname, 'g', 'H', 'triGroup', 'shapes', 'domain', 'contour');
            fprintf('[MeshBuilder] Maillage sauvegardé : %s\n', fname);
        end

        function show(obj)
            % SHOW  Affiche le maillage coloré par groupes + contour extérieur.
            %
            % - Patches colorés par étiquette `triGroup`
            % - Contour extérieur superposé (rouge)
            % - Légende via `obj.group_names`
            %
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
