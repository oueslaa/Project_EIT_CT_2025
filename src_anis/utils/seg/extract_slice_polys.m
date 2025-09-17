function [domain, shapes, ext_smooth, info] = extract_slice_polys(patient_id, z_slice, opts)
% EXTRACT_SLICE_POLYS  Extrait le contour externe + polyshapes d'organes pour une slice CT.
%
% But
% ----
% À partir du CT et des segmentations 3D d'un patient, produire :
%  - ext_smooth : contour extérieur (mm) lissé/échantillonné
%  - shapes     : polyshapes par grands groupes d'organes (Heart, Lung, ...)
%  - domain     : "soft tissue" = extérieur - (union des organes)
%
% Entrées
% -------
% patient_id : ex. 's0011'  (doit exister dans Data_set/)
% z_slice    : indice de coupe (1..D, D = taille Z du CT)
% opts       : (optionnel) struct de réglages :
%    .ct_filename     (def 'ct.nii.gz')
%    .seg_dirname     (def 'segmentations')
%    .ct_thresh_hu    (def -400) seuillage pour isoler le corps
%    .min_area_mm2    (def 300)  surface min d'un blob à garder
%    .smooth_K        (def 50)   lissage FFT du contour externe
%    .smooth_step     (def 5)    sous-échantillonnage externe
%    .smooth_K_org    (def 30)   lissage FFT des organes
%    .smooth_step_org (def 2)    sous-échantillonnage organes
%
% Sorties
% -------
% domain     : polyshape (extérieur – organes, nettoyé)
% shapes     : struct de polyshape par groupe (Heart, Lung, Trachea, Bone, ...)
% ext_smooth : [N x 2] contour extérieur lissé (mm, fermé)
% info       : niftiinfo du CT
%
% Dépendances
% -----------
% - read_nifti3D, fftSmooth, rmholes
%
% Pièges & bonnes pratiques
% -------------------------
% - La matrice d'affine (voxel->mm) est lue depuis info.Transform.T quand
%   disponible, sinon fallback diag(PixelDimensions) (pas de translation).
% - Les poumons : on garde les 2 plus gros blobs (gauche/droite) par défaut.
% - Le "priority peeling" évite des recouvrements (Heart > Lung > Bone > ...).
% - Ajuste ct_thresh_hu selon le CT (ex. [-500,-300]) si le corps sort mal.
%
% Exemple
% -------
% [domain,shapes,ext,info] = extract_slice_polys('s0011', 310);
%
% ©2025
%
% ---------- Options par défaut ----------
if nargin < 3, opts = struct(); end
def = struct('ct_filename','ct.nii.gz','seg_dirname','segmentations', ...
             'ct_thresh_hu',-400,'min_area_mm2',300, ...
             'smooth_K',50,'smooth_step',5, ...
             'smooth_K_org',30,'smooth_step_org',2);
fn = fieldnames(def);
for i=1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = def.(fn{i}); end
end

% ---------- Localisation des fichiers ----------
data_dir = fullfile('Data_set', patient_id);
ct_file  = fullfile(data_dir, opts.ct_filename);
seg_dir  = fullfile(data_dir, opts.seg_dirname);
assert(isfolder(data_dir), 'Patient folder not found: %s', data_dir);
assert(exist(ct_file,'file')==2, 'CT file not found: %s', ct_file);
assert(isfolder(seg_dir), 'Seg folder not found: %s', seg_dir);

% ---------- Lecture CT ----------
[CT, info] = read_nifti3D(ct_file);
[H,W,D] = size(CT);
assert(z_slice>=1 && z_slice<=D, 'z_slice out of bounds (1..%d)', D);

% ---------- Affine voxel->mm ----------
T = local_get_affine(info);

% ---------- 1) Contour extérieur ----------
I  = CT(:,:,z_slice);
bw = I > opts.ct_thresh_hu;

% Morphologie robuste (taille strel basée sur la résolution in-plane)
px = local_pixdim(info);
se_radius = max(1, round(8 / max(px(1),eps)));
bw = imclose(bw, strel('disk', se_radius));
bw = imfill(bw, 'holes');
bw = bwareaopen(bw, 10000);

CC = bwconncomp(bw,8);
assert(CC.NumObjects>0, 'No body component detected on slice %d.', z_slice);
np  = cellfun(@numel, CC.PixelIdxList);
[~,mi] = max(np);
body = false(size(bw)); body(CC.PixelIdxList{mi}) = true;

B = bwboundaries(body,'noholes');
ext_pix = B{1}(:,[2,1]); % [x_pix,y_pix]

% pixel->mm pour la slice z
Npt   = size(ext_pix,1);
XYZ1  = [ ext_pix(:,1)'; ext_pix(:,2)'; repmat(z_slice,1,Npt); ones(1,Npt) ];
XYZmm = (T * XYZ1)';                      % Nx4
ext_mm = XYZmm(:,1:2);

% lissage et fermeture
ext_smooth = fftSmooth(ext_mm, opts.smooth_K, opts.smooth_step);
if ~isequal(ext_smooth(1,:), ext_smooth(end,:))
    ext_smooth(end+1,:) = ext_smooth(1,:);
end

% ---------- 2) Chargement segmentations ----------
S = dir(fullfile(seg_dir,'*.nii*'));
volumes = struct();
for k=1:numel(S)
    nm   = S(k).name;
    base = regexprep(nm,'\.nii(\.gz)?$','');
    fld  = matlab.lang.makeValidName(base);
    V3   = read_nifti3D(fullfile(seg_dir,nm));
    if size(V3,3) < z_slice, continue; end
    volumes.(fld) = squeeze(V3(:,:,z_slice)) > 0;
end
fnv = fieldnames(volumes);

% ---------- 3) Mapping nom -> groupe ----------
organ_natures = struct( ...
  'Heart',                 {{'heart','atrial_appendage'}}, ...
  'Lung',                  {{'lung','pulmonary_lobe','upper_lobe','middle_lobe','lower_lobe'}}, ...
  'Trachea',               {{'trachea'}}, ...
  'Bone',                  {{'rib','vertebrae','vertebrae_','clavicle','clavicula','scapula','sternum','skull','humerus','femur','sacrum','hip'}}, ...
  'Cartilage',             {{'costal_cartilage','costal_cartilages','cartilage'}}, ...
  'Blood',                 {{'aorta','vena_cava','venaecava','inferior_vena_cava','superior_vena_cava','brachiocephalic','subclavian','carotid','common_carotid','iliac_artery','iliac_vena','pulmonary_vein','portal_vein','splenic_vein'}}, ...
  'LiverSpleenPancreas',   {{'liver','spleen','pancreas'}}, ...
  'Kidney',                {{'kidney','kidney_left','kidney_right','kidney_cyst'}}, ...
  'Other',                 {{'esophagus','stomach','small_bowel','duodenum','colon','gallbladder','adrenal_gland','thyroid_gland','spinal_cord','brain','prostate','autochthon','iliopsoas','gluteus','hip_left','hip_right'}} ...
);

organ_group = containers.Map();
for i=1:numel(fnv)
    nm = fnv{i}; low = lower(nm); grp = 'Other';
    gnames = fieldnames(organ_natures);
    for gi = 1:numel(gnames)
        toks = organ_natures.(gnames{gi});
        hit = false;
        for tj = 1:numel(toks)
            if contains(low, toks{tj})
                grp = gnames{gi}; hit = true; break;
            end
        end
        if hit, break; end
    end
    organ_group(nm) = grp;
end

% ---------- 4) Fusion/lissage par groupe (en mm) ----------
groups = {'Heart','Lung','Trachea','Bone','Cartilage','Blood','LiverSpleenPancreas','Kidney','Other'};
shapes = struct(); for ig=1:numel(groups), shapes.(groups{ig}) = polyshape(); end

% poumons : on conserve les 2 plus grands blobs
for ig = 1:numel(groups)
    gname = groups{ig}; %#ok<NBRAK>  % <- MATLAB autorise {}, on garde la forme standard ci-dessous
end
for ig = 1:numel(groups)
    gname = groups{ig};
    if strcmp(gname,'Lung')
        mask_lung = false(H,W);
        for i=1:numel(fnv)
            nm = fnv{i};
            if strcmp(organ_group(nm),'Lung')
                mask_lung = mask_lung | volumes.(nm);
            end
        end
        if any(mask_lung(:))
            CC = bwconncomp(mask_lung,8);
            npix = cellfun(@numel, CC.PixelIdxList);
            kkeep = min(2, numel(npix));
            [~,idx2] = maxk(npix, kkeep);
            for j = 1:numel(idx2)
                mask_j = false(H,W); mask_j(CC.PixelIdxList{idx2(j)}) = true;
                B = bwboundaries(mask_j,'noholes');
                for b=1:numel(B)
                    pts_pix = B{b}(:,[2,1]); Np = size(pts_pix,1);
                    if Np < 3, continue; end
                    pts_vox = [pts_pix, repmat(z_slice, Np,1)];
                    XYZ1 = [pts_vox, ones(Np,1)]';
                    pts_mm = (T * XYZ1)'; pts_mm = pts_mm(:,1:2);
                    pts_mm = fftSmooth(pts_mm, opts.smooth_K_org, opts.smooth_step_org);
                    if size(pts_mm,1) < 3, continue; end
                    P = polyshape(pts_mm(:,1), pts_mm(:,2), 'Simplify', true);
                    shapes.Lung = union(shapes.Lung, P);
                end
            end
        end
        continue
    end

    % autres groupes
    for i=1:numel(fnv)
        nm = fnv{i};
        if ~strcmp(organ_group(nm), gname), continue; end
        mask = volumes.(nm);
        if ~any(mask(:)), continue; end
        B = bwboundaries(mask,'holes');
        for b=1:numel(B)
            pts_pix = B{b}(:,[2,1]); Np = size(pts_pix,1);
            if Np < 3, continue; end
            pts_vox = [pts_pix, repmat(z_slice,Np,1)];
            XYZ1 = [pts_vox, ones(Np,1)]';
            pts_mm = (T * XYZ1)'; pts_mm = pts_mm(:,1:2);
            pts_mm = fftSmooth(pts_mm, opts.smooth_K_org, opts.smooth_step_org);
            if size(pts_mm,1) < 3, continue; end
            if abs(polyarea(pts_mm(:,1), pts_mm(:,2))) < opts.min_area_mm2, continue; end
            P = polyshape(pts_mm(:,1), pts_mm(:,2), 'Simplify', true);
            shapes.(gname) = union(shapes.(gname), P);
        end
    end
end

% ---------- 5) Priorité d'intersection (cohérence anatomique simple) ----------
priority = {'Heart','Lung','Bone','Cartilage','Blood','LiverSpleenPancreas','Kidney','Other'};
for i = 2:numel(priority)
    for j = 1:i-1
        shapes.(priority{i}) = subtract(shapes.(priority{i}), shapes.(priority{j}));
        shapes.(priority{i}) = rmholes(shapes.(priority{i}));
    end
end

% ---------- 6) Soft tissue (extérieur - organes) ----------
domain = polyshape(ext_smooth(:,1), ext_smooth(:,2));
for ig = 1:numel(priority)
    domain = subtract(domain, shapes.(priority{ig}));
    domain = simplify(domain, 'KeepCollinear', true);
    domain = rmholes(domain);
end
domain = simplify(domain,'KeepCollinear',true);
domain = rmholes(domain);

end

% ====== Helpers locaux (privés à ce fichier) =============================

function T = local_get_affine(info)
% Essaie info.Transform.T, sinon fallback basique diag(PixelDimensions)
    if isfield(info,'Transform') && isfield(info.Transform,'T')
        T = info.Transform.T;
    else
        % Fallback minimal : diag des pixel sizes, sans translation
        pd = [1 1 1];
        if isfield(info,'PixelDimensions') && numel(info.PixelDimensions)>=3
            pd = double(info.PixelDimensions(1:3));
        end
        T = eye(4);
        T(1,1)=pd(1); T(2,2)=pd(2); T(3,3)=pd(3);
    end
end

function pd = local_pixdim(info)
% Récupère la résolution in-plane [dx dy] (mm)
    if isfield(info,'PixelDimensions') && numel(info.PixelDimensions)>=2
        pd = double(info.PixelDimensions(1:2));
    else
        pd = [1 1];
    end
end
