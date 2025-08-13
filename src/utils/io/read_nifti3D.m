function [V,info] = read_nifti3D(fp)
% Lit un .nii ou .nii.gz proprement sans collision de fichiers temporaires
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