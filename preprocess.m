% preprocess.m
% Pipeline de traitement des segmentations NIfTI, fusion de masques et export.

%% 1. PARAMÈTRES UTILISATEUR
case_id        = 's0011';
z_min          = 201;
z_max          = 210;
seg_dir        = fullfile('Data_set', case_id, 'segmentations');
output_seg_dir = 'segmentation_slices';
merged_dir     = 'merged_masks';

if ~exist(output_seg_dir, 'dir'), mkdir(output_seg_dir); end
if ~exist(merged_dir, 'dir'),    mkdir(merged_dir);    end

addpath('src');  % Assure que Create_skin_mask_bis et Create_organ_mask sont accessibles

%% 2. Dictionnaire pour noms courts
short_names = containers.Map( ...
    {'duodenum','pancreas','stomach','autochthon_left','autochthon_right', ...
     'vertebrae_L1','vertebrae_L2','costal_cartilages','rib_left_10', ...
     'rib_left_11','rib_left_12','rib_right_11','rib_right_12','liver', ...
     'aorta','gallbladder','inferior_vena_cava','small_bowel','colon', ...
     'kidney_left','kidney_right','iliopsoas_left','iliopsoas_right', ...
     'spinal_cord','soft_tissue'}, ...
    {'duod','panc','stom','autoL','autoR','vertL1','vertL2','cost', ...
     'ribL10','ribL11','ribL12','ribR11','ribR12','live','aort', ...
     'gall','ivc','sbow','colo','kidL','kidR','ipsL','ipsR','cord','soft'} );
short = @(name) ( ...
    short_names.isKey(name) * string(short_names(name)) + ...
   ~short_names.isKey(name) * string(name(1:min(4,end))) );

%% 3. LISTER LES FICHIERS ET CONSTRUIRE ORGANS (map key->1x1 cell)
d = dir(fullfile(seg_dir,'*.nii.gz'));
ORGANS = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(d)
    fname = d(i).name;
    key   = erase(fname, '.nii.gz');
    ORGANS(key) = { fname };   % stocke dans une cellule
end

%% 4. DÉTECTER LES ORGANES PRÉSENTS SUR LA PLAGE Z
organs_presnt = organs_present_in_crop(seg_dir, z_min, z_max, ORGANS);

%% 5. CHARGEMENT DES MASQUES BINAIRES
masks = struct();
for i = 1:numel(organs_presnt)
    organ = organs_presnt{i};
    cellVal = ORGANS(organ);
    fname   = cellVal{1};  % extraction de la string
    fullp   = fullfile(seg_dir, fname);
    info    = niftiinfo(fullp);
    vol     = niftiread(info);
    crop3d  = vol(:,:, z_min:z_max);
    masks.(organ).data = crop3d > 0;
    masks.(organ).info = info;
end

%% 6. AJOUT DU MASQUE SOFT_TISSUE
mask_cells = Get_soft_tissue_mask(case_id, organs_presnt, ORGANS);
soft3d     = mask_cells{1};                    % soft_tissue brut
soft_slice = soft3d(:,:, z_min:z_max) > 0;      % même tranche Z
ref_info   = masks.(organs_presnt{1}).info;
masks.soft_tissue.data = soft_slice;
masks.soft_tissue.info = ref_info;

%% 7. DÉTECTION DES COUPLES QUI SE TOUCHENT
keys   = fieldnames(masks);
nSlice = size(masks.(keys{1}).data, 3);
pairs  = strings(0,2);
for z = 1:nSlice
    for i = 1:numel(keys)
        for j = i+1:numel(keys)
            m1 = masks.(keys{i}).data(:,:,z);
            m2 = masks.(keys{j}).data(:,:,z);
            if any(m1 & m2, 'all')
                pair = sort([ string(keys{i}), string(keys{j}) ]);
                pairs(end+1,:) = pair; %#ok<SAGROW>
            end
        end
    end
end
if ~isempty(pairs)
    pairs = unique(pairs, 'rows');
end

%% 8. FUSION DES MASQUES TOUCHANTS
used = strings(0,1);
merged_masks = struct();
merged_map   = struct();
for k = 1:size(pairs,1)
    a = pairs(k,1); b = pairs(k,2);
    if any(used==a) || any(used==b), continue; end
    union3d = masks.(char(a)).data | masks.(char(b)).data;
    key_short = sprintf('%s_%s', short(char(a)), short(char(b)));
    merged_masks.(key_short).data = union3d;
    merged_masks.(key_short).info = masks.(char(a)).info;
    merged_map.(key_short) = { char(a), char(b) };
    used = [used; a; b];
end
for i = 1:numel(keys)
    kname = keys{i};
    if ~any(used==string(kname))
        ks = short(kname);
        merged_masks.(char(ks)) = masks.(kname);
        merged_map.(char(ks))     = { kname };
    end
end

%% 9. ATTRIBUTION DES COULEURS (colormap 'lines')
mkeys = fieldnames(merged_masks);
nKeys = numel(mkeys);
cmap  = lines(nKeys);
group_colors = containers.Map('KeyType','char','ValueType','any');
for i = 1:nKeys
    rgba = [ cmap(i,:), 1 ];  % alpha=1
    group_colors(mkeys{i}) = rgba;
end

%% 10. SAUVEGARDE ET AFFICHAGE PAR SLICE
h = size(merged_masks.(mkeys{1}).data, 1);
w = size(merged_masks.(mkeys{1}).data, 2);

for zIdx = z_min:z_max
    slice_loc = zIdx - z_min + 1;
    seg_img = zeros(h, w, 4);
    for i = 1:nKeys
        key    = mkeys{i};
        vol3d  = merged_masks.(key).data;
        mask2d = vol3d(:,:,slice_loc);
        color  = group_colors(key);
        for c = 1:4
            layer = seg_img(:,:,c);
            layer(mask2d) = color(c);
            seg_img(:,:,c) = layer;
        end
    end
    RGB   = seg_img(:,:,1:3);
    alpha = seg_img(:,:,4);
    outfn = fullfile(output_seg_dir, sprintf('segmented_slice_%d.png', zIdx));
    imwrite(RGB, outfn, 'PNG', 'Alpha', alpha);
end
disp(['✅ Images segmentées dans "', output_seg_dir, ...
      '" pour les slices ', num2str(z_min), ' à ', num2str(z_max)]);

%% 11. SAUVEGARDE DU CODE COULEUR (.mat)
group_names = mkeys;
rgb_colors  = zeros(nKeys,3);
for i = 1:nKeys
    tmp = group_colors(mkeys{i});
    rgb_colors(i,:) = tmp(1:3);
end
save('group_colors.mat', 'group_names', 'rgb_colors');
disp('✅ Couleurs sauvegardées dans "group_colors.mat"');

%% 12. SAUVEGARDE DES MASQUES FUSIONNÉS + MAPPING JSON
for i = 1:nKeys
    key  = mkeys{i};
    data = merged_masks.(key).data;
    info = merged_masks.(key).info;
    outfn = fullfile(merged_dir, sprintf('%s_mask_z%d_%d.nii.gz', ...
                   key, z_min, z_max));
    niftiwrite(uint8(data), outfn, info, 'Compressed', true);
end
json_txt = jsonencode(merged_map);
fid = fopen(fullfile(merged_dir,'merged_mapping.json'), 'w');
fprintf(fid, '%s', json_txt);
fclose(fid);
disp(['✅ Masques fusionnés et mapping dans "', merged_dir, '"']);

%% Fonctions auxiliaires

function present = organs_present_in_crop(seg_dir, z_min, z_max, ORGANS)
    keys = ORGANS.keys;
    present = {};
    for i = 1:numel(keys)
        organ = keys{i};
        parts = ORGANS(organ);  % cell array
        found = false;
        for j = 1:numel(parts)
            segp = fullfile(seg_dir, parts{j});
            if ~exist(segp,'file'), continue; end
            info = niftiinfo(segp);
            vol  = niftiread(info);
            if any(vol(:,:,z_min:z_max),'all')
                found = true; break;
            end
        end
        if found
            present{end+1} = organ; %#ok<AGROW>
        end
    end
end

function mask3d = Get_soft_tissue_mask(case_id, organs_presnt, ORGANS)
    selected = {};
    for i = 1:numel(organs_presnt)
        parts = ORGANS(organs_presnt{i});
        selected = [selected, parts]; %#ok<AGROW>
    end
    [soft_tissue, outside, ~] = Create_skin_mask_bis(case_id, selected);
    mask3d = {soft_tissue, outside};
    for i = 1:numel(organs_presnt)
        organ_parts = ORGANS(organs_presnt{i});
        organ_mask  = Create_organ_mask(case_id, organ_parts);
        mask3d{end+1} = organ_mask; %#ok<AGROW>
    end
end
