%% Example real data
%% First you need to add the pathes that contains all the needed functions

clear all
path = '/Users/hao/Dropbox/DMRI_code/';

addpath(path);
addpath(strcat(path,'ADMM'));
addpath(strcat(path,'construction_functions'));
addpath(strcat(path,'toolbox_wavelet_meshes'));
addpath(strcat(path,'toolbox_wavelet_meshes/toolbox/'));
addpath(strcat(path,'MEALPix/'));
addpath(strcat(path,'help_functions/'));
addpath(strcat(path,'NIfTI/'));

addpath('/Users/hao/Dropbox/stats_project/FOD_codes_simulation/matrix');
path_save='/Users/hao/Dropbox/stats_project/FOD_codes_simulation/Real_data/';
subject_id = 'mgh_1022';
path_folder = strcat(path_save, subject_id, '/diff/preproc/');

addpath(path_folder);
addpath(strcat(path_folder,'mri/'));

%% load nii data (raw)
bvec_raw = load(strcat(path_folder,'bvecs_fsl_moco_norm.txt')); 
bvec_raw = bvec_raw';
bval_raw = load(strcat(path_folder,'bvals.txt'));
bval_raw = bval_raw';
nii_data = load_untouch_nii(strcat(path_folder,'mri/diff_preproc.nii.gz'));


b = 5000;
b_sample = find(abs(bval_raw-b)==0);
b0_sample = find(abs(bval_raw)==0);
b_idx = sort([b_sample, b0_sample]);

img_data_all = nii_data.img(:,:,:,b_idx);

nii_data_new = nii_data;
nii_data_new.img = img_data_all;
nii_data_new.fileprefix = strcat(path_folder,'nii_b',num2str(b/1000));
nii_data_new.hdr.dime.dim(5) = length(b_idx);
bvec = bvec_raw(:,b_idx);
bval = bval_raw(:,b_idx);

save_path = strcat(path_folder,'nii_b',num2str(b/1000),'/'); 
mkdir(save_path);
save_untouch_nii(nii_data_new, strcat(save_path,'nii_data_b',num2str(b/1000),'.nii.gz'));
save(strcat(save_path,'bvecs.txt'),'bvec','-ascii');
save(strcat(save_path,'bvals.txt'),'bval','-ascii');

