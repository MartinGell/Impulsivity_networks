function parcellate(img, mask, outputdir)
% img = path/name of nii image
% mask = path to mask nii
%cd('/home/mgell/Matlab_toolboxes/Juspace_v1/PETatlas');

%outputdir = '/home/mgell/Work/FC/parcellated_PET_extended/';
%mask = '/home/mgell/Downloads/BN_Atlas_246_3mm.nii.gz'
%mask = '/home/mgell/Work/FC/nets/imp_all_net_3x3x3.nii';
%mask = '/home/mgell/Work/FC/nets/power_nodes_3x3x3.nii';
%mask =  '/home/mgell/Work/FC/nets/additional_regions/Extended_net_overlap_RL_3x3x3.nii';

%addpath('/home/mgell/Work/FC/PhD/');
addpath('/home/mgell/Matlab_toolboxes/spm12');
addpath('/home/mgell/Matlab_toolboxes/NIFTI_toolbox-master');
addpath(genpath('/home/mgell/Matlab_toolboxes/DPABI_V5.0_201001/'));
addpath(genpath('/home/mgell/Work/FC/code/'));


nii = spm_vol(img);
nii_data = spm_read_vols(nii);

X = spm_vol(mask);
assert(all(all(nii.mat == X.mat)), 'Different orientation');
assert(all(all(nii.dim == X.dim)), 'Different size');

img = strsplit(img, '.');
extract_roi_tc(nii_data, mask, 2, [outputdir img{1} '_iso2_tri_CAT.txt']);
%extract_roi_tc(nii_data, mask, 2, [outputdir 'Power_' img '.txt']);
%extract_roi_tc(nii_data, mask, 2, outputdir);

end