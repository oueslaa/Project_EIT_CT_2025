format compact;
addpath src_aatae/

file_subject = 'Data_set/s0011/ct.nii.gz';
file_organ = 'Data_set/s0011/segmentations/liver.nii.gz';

% visualise_nii(file_subject);  
% visualise_cube(file_organ);  

wrap = visualise_alpha(file_organ);  

wrap.volume
wrap.surfaceArea
