%% OOEIT - Setup complet (Forward + Inverse) aligné avec [Howard et al., 2025]
clear; close all; clc;

%% --- Paramètres généraux ---
z_slice    = 301;
meshFile   = fullfile('MeshData', sprintf('mesh_slice%03d.mat', z_slice));
Ne         = 16;          % Nombre d'électrodes
el_width   = 33;          % Largeur électrode en mm (3.3 cm)
I_amp      = 0.35e-3;     % Courant injecté (A)
z_contact  = 0.05;        % Ω·m² (contact impedance)

% === Conductivités par groupe (proches du papier, Table 1) ===
% Paper : Heart 0.5, SoftTissue 0.3, Lung 0.15, Trachea 0.15, Bone 0.05
condMap = containers.Map( ...
    {'SoftTissue','Heart','Lung','Trachea','Bone','Other'}, ...
    [0.3,          0.5,     0.15,   0.15,     0.05,  0.3] ...
);

%% --- Charger mesh et infos ---
load(meshFile, 'g', 'H', 'triGroup', 'domain');
Mtri = size(H,1);

%% --- 1) Conductivité par triangle (sans perturbation aléatoire pour reproductibilité) ---
sigma_tri = zeros(Mtri,1);
types = keys(condMap);
for t = 1:Mtri
    grpIdx = triGroup(t);
    if grpIdx > 0 && grpIdx <= numel(types)
        sigma_tri(t) = condMap(types{grpIdx});
    else
        sigma_tri(t) = condMap('SoftTissue');
    end
end

%% --- 2) Placement des électrodes (arêtes de bord) ---
extVerts = domain.Vertices;
[E, el_centers] = buildElectrodesFromContour(g, extVerts, Ne, el_width, H);


g = g / 1000;              % mm → m
extVerts = extVerts / 1000;
el_centers = el_centers / 1000; % (si tu les utilises pour visualisation)

%% --- 3) Création du mesh OOEIT ---
mesh = ForwardMesh1st(g, H, E);

%% --- 4) Solveur EIT (CEM, courant pattern trigonométrique) ---
solver = EITFEM(mesh);
solver.zeta = z_contact*ones(Ne,1);   % Impédance de contact
solver.mode = 'current';

% Pattern trigonométrique (voir papier)
injPattern = buildTrigPattern(Ne, I_amp);
solver.Iel = injPattern(:);  % Injection pour le forward

%% --- 5) Résolution Forward (simulation des mesures) ---
sigma_used = solver.PreProcessSigma(sigma_tri);
Imeas_clean = solver.SolveForwardVec(sigma_tri);

% === Ajout du bruit gaussien 0.1% du signal max (cf. papier, Sec 2.1.2) ===
noise = 0.001 * max(abs(Imeas_clean)) * randn(size(Imeas_clean));
Imeas = Imeas_clean + noise;

%% --- 6) Affichages Forward ---
% Figure 1 : Mesh + électrodes
figure('Name','Mesh + électrodes','Color','w');
triplot(H, g(:,1), g(:,2), 'Color', [0.6 0.6 0.6]); hold on;
plot([extVerts(:,1); extVerts(1,1)], [extVerts(:,2); extVerts(1,2)], 'r-', 'LineWidth', 1.5);
scatter(el_centers(:,1), el_centers(:,2), 50, 'filled', 'r');
for k = 1:Ne
    text(el_centers(k,1), el_centers(k,2), sprintf('%d',k), ...
        'FontSize',10, 'FontWeight','bold', 'Color','b', ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
axis equal;
title('Mesh + électrodes (numérotées)');

% Figure 2 : σ utilisée
figure('Name','\sigma utilisée','Color','w');
patch('Faces',H,'Vertices',g,'FaceVertexCData',sigma_used,...
      'FaceColor','flat','EdgeColor','none');
axis equal tight; colormap(jet); colorbar;
title('\sigma utilisée par le solveur [S/m]');

% Figure 3 : Signaux simulés (avec bruit)
figure('Name','Signaux simulés','Color','w');
plot(Imeas,'-o','LineWidth',1.5);
xlabel('Index mesure'); ylabel('Tension [V]');
title('Signaux simulés (pattern trigonométrique, bruit 0.1%)');
grid on;

%% --- 7) Problème inverse ---
TVPrior = PriorTotalVariation(g, H, 3); % (lambda à ajuster si besoin)
ofuns = {solver; TVPrior};
invSolver = SolverGN(ofuns);
invSolver.maxIterInLine = 15;
invSolver.plotter = Plotter(g, H);

% On donne les mesures simulées au solver (mode = 'current' → Uel = potentiels)
solver.Uel = Imeas;

% Reconstruction avec un guess homogène
sig_init = ones(size(g,1),1) * condMap('SoftTissue');
sigma_rec = invSolver.Solve(sig_init);

%% --- 8) Affichage reconstruction ---
figure('Name','Reconstruction finale','Color','w');
patch('Faces',H,'Vertices',g,'FaceVertexCData',sigma_rec,...
      'FaceColor','flat','EdgeColor','none');
axis equal tight; colormap(jet); colorbar;
title('Conductivité reconstruite [S/m]');

%% --- Fonctions utilitaires (inchangées, voir ton code plus haut) ---
% buildElectrodesFromContour
% projectPointOnContour
% closestPointOnSegment
% buildTrigPattern

% [Fonctions utilitaires collées ici si besoin]


%% --- Fonctions utilitaires ---
function [E, el_centers] = buildElectrodesFromContour(g, contour, Ne, el_width, H)
    d = sqrt(sum(diff(contour).^2,2));
    arcLength = [0; cumsum(d)];
    perim = arcLength(end);
    spacing = perim / Ne;

    edges_all = unique(sort([H(:,[1 2]); H(:,[2 3]); H(:,[3 1])],2),'rows');
    edge_count = zeros(size(edges_all,1),1);
    for i = 1:size(edges_all,1)
        e = edges_all(i,:);
        count = sum(ismember(sort([H(:,[1 2]); H(:,[2 3]); H(:,[3 1])],2), e, 'rows'));
        edge_count(i) = count;
    end
    border_edges = edges_all(edge_count==1,:);

    proj_pos = zeros(size(g,1),1);
    for n = 1:size(g,1)
        [~, proj_pos(n)] = projectPointOnContour(g(n,:), contour, arcLength);
    end

    E = cell(Ne,1);
    el_centers = zeros(Ne,2);

    for k = 1:Ne
        center_pos = (k-1)*spacing + spacing/2;
        idx = find(arcLength <= center_pos, 1, 'last');
        ratio = (center_pos - arcLength(idx)) / (arcLength(idx+1)-arcLength(idx));
        el_centers(k,:) = contour(idx,:) + ratio*(contour(idx+1,:) - contour(idx,:));

        arc_start = mod(center_pos - el_width/2, perim);
        arc_end   = mod(center_pos + el_width/2, perim);

        if arc_start < arc_end
            in_arc = @(pos) pos >= arc_start & pos <= arc_end;
        else
            in_arc = @(pos) pos >= arc_start | pos <= arc_end;
        end

        in_elec = arrayfun(@(n) in_arc(proj_pos(n)), 1:size(g,1));
        mask_edges = all(in_elec(border_edges), 2);
        E{k} = border_edges(mask_edges,:);
    end
end

function [dist, pos] = projectPointOnContour(pt, contour, arcLength)
    minDist = inf; posOnArc = 0;
    for i = 1:size(contour,1)-1
        seg = [contour(i,:); contour(i+1,:)];
        [p_proj, t] = closestPointOnSegment(pt, seg);
        dist_tmp = norm(pt - p_proj);
        if dist_tmp < minDist
            minDist = dist_tmp;
            posOnArc = arcLength(i) + t*norm(diff(seg));
        end
    end
    dist = minDist;
    pos = posOnArc;
end

function [p_proj, t] = closestPointOnSegment(p, seg)
    v = seg(2,:) - seg(1,:);
    w = p - seg(1,:);
    t = dot(w,v) / dot(v,v);
    t = max(0,min(1,t));
    p_proj = seg(1,:) + t*v;
end

function injPattern = buildTrigPattern(Ne, I_amp)
    injPattern = zeros(Ne,Ne);
    for k = 1:Ne
        phase = 2*pi*(k-1)/Ne;
        for n = 1:Ne
            injPattern(n,k) = I_amp * sin( 2*pi*(n-1)/Ne + phase );
        end
    end
end
