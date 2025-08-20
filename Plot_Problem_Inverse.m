% ===== Plot_Problem_Inverse.m =====
clear; close all; clc; addpath('src', genpath('src'));

patient_id = 's0011'; z_slice = 301;

rootOut = fullfile('Outputs', patient_id, sprintf('slice_%03d', z_slice));
recFile = fullfile(rootOut, 'reconstruction_inverse.mat');
packFile= fullfile(rootOut, 'eit_pack.mat');
vizFile = fullfile(rootOut, 'viz_config.mat');
plotDir = fullfile(rootOut, 'plots_inverse'); if ~exist(plotDir,'dir'), mkdir(plotDir); end
assert(isfile(recFile), 'Reconstruction manquante: %s', recFile);

R = load(recFile);
S = []; if isfile(packFile), S = load(packFile); end
viz = []; if isfile(vizFile), tmp = load(vizFile); viz = tmp.viz_config; end

if ~isempty(S), g=S.g; H=S.H; else, g=R.g; H=R.H; end
tri0   = node2tri_avg_local(H, R.sigma0);
triRec = R.sigma_rec_tri; triGT=[];
if isfield(R,'sigma_true_tri') && ~isempty(R.sigma_true_tri), triGT=R.sigma_true_tri(:);
elseif ~isempty(S) && isfield(S,'sigma_tri'), triGT = S.sigma_tri(:); end

climUsed = []; if ~isempty(viz) && isfield(viz,'sigma_clim'), climUsed = viz.sigma_clim; end

fh = figure('Color','w','Name','Triptyque','NumberTitle','off');
tl = tiledlayout(fh,1,3,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
patch('Faces',H,'Vertices',g,'FaceVertexCData',tri0,'FaceColor','flat','EdgeColor','none');
axis equal tight; colormap(turbo); if ~isempty(climUsed), caxis(climUsed); end
title('Init (\sigma_0)'); cb=colorbar; ylabel(cb,'S/m');

nexttile; hold on;
if ~isempty(triGT)
  patch('Faces',H,'Vertices',g,'FaceVertexCData',triGT,'FaceColor','flat','EdgeColor','none');
  if ~isempty(climUsed), caxis(climUsed); end
  title('Cible (\sigma_{true})'); cb=colorbar; ylabel(cb,'S/m');
else
  text(0.5,0.5,'Ground truth indisponible','Units','normalized','HorizontalAlignment','center'); axis off;
end

nexttile; hold on;
patch('Faces',H,'Vertices',g,'FaceVertexCData',triRec,'FaceColor','flat','EdgeColor',[.6 .6 .6],'EdgeAlpha',0.25);
axis equal tight; colormap(turbo); if ~isempty(climUsed), caxis(climUsed); end
title(sprintf('Reco — misfit=%.4g', R.misfit)); cb=colorbar; ylabel(cb,'S/m');

exportgraphics(fh, fullfile(plotDir, sprintf('triptyque_slice%03d.png',z_slice)), 'Resolution', 300);

function tri = node2tri_avg_local(H, nod)
    nod = nod(:);
    tri = (nod(H(:,1)) + nod(H(:,2)) + nod(H(:,3))) / 3;
    tri = tri(:);
end
