function [V,info] = read_nifti3D(fp)
% READ_NIFTI3D  Lit un .nii ou .nii.gz proprement (sans collision de fichiers temp).
%
% But
% ----
% Charger un volume NIfTI 3D (ex: CT, masque) et renvoyer à la fois les
% voxels et le header NIfTI (niftiinfo).
%
% Entrée
% ------
% fp : chemin fichier .nii ou .nii.gz
%
% Sorties
% -------
% V    : volume 3D (type et taille selon le contenu)
% info : struct NIfTI (niftiinfo)
%
% Remarques
% ---------
% - Pour .nii.gz, on décompresse dans un dossier temporaire unique puis on
%   nettoie à la fin via onCleanup.
% - Ne modifie pas le fichier d'entrée.
%
% Exemple
% -------
% [CT,info] = read_nifti3D('Data_set/s0011/ct.nii.gz');
%
    if endsWith(fp,'.nii.gz')
        tmpdir = tempname; mkdir(tmpdir);
        cleaner = onCleanup(@() rmdir(tmpdir,'s')); %#ok<NASGU>
        gunzip(fp, tmpdir);
        [~,b,~] = fileparts(fp);
        tmp = fullfile(tmpdir,b);
        info = niftiinfo(tmp);
        V    = niftiread(tmp);
    else
        info = niftiinfo(fp);
        V    = niftiread(fp);
    end
end
