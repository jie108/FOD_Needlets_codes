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
subject_id = 'S110933';
path_folder = strcat(path_save, subject_id);

addpath(path_folder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load .nii data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bvec_raw = load('027-S-2245_060611_S110933.bvec'); 
bval_raw = load('027-S-2245_060611_S110933.bval');
% nii_data = load_nii(strcat(path_folder,'027-S-2245_060611_S110933_data.nii.gz'));
% nii_data_header = nii_data.hdr;

%% load fsl processed data
nii_FA = load_nii('dti_FA.nii.gz');
nii_MO = load_nii('dti_MO.nii.gz');
nii_MD = load_nii('dti_MD.nii.gz');
nii_S0 = load_nii('dti_S0.nii.gz');

nii_V1 = load_nii('dti_V1.nii.gz');
nii_V2 = load_nii('dti_V2.nii.gz');
nii_V3 = load_nii('dti_V3.nii.gz');

nii_L1 = load_nii('dti_L1.nii.gz');
nii_L2 = load_nii('dti_L2.nii.gz');
nii_L3 = load_nii('dti_L3.nii.gz');

JHU_wm_label = load_nii('nat_JHU_wm_label.nii.gz');
ATR_R = load_nii('nat_Anterior_thalamic_radiation_R.nii.gz');
Cigulum = load_nii('nat_Cingulum_(cingulate_gyrus)_R.nii.gz');
CT_R = load_nii('nat_Corticospinal_tract_R.nii.gz');
Fm_R = load_nii('nat_Forceps_minor.nii.gz');
IFOF_R = load_nii('nat_Inferior_fronto-occipital_fasciculus_R.nii.gz');
SLF_R = load_nii('nat_Superior_longitudinal_fasciculus_R.nii.gz');

ROI = load_nii('ROI.nii.gz');

%% transfer to matrix
% img_data_raw = nii_data.img;
% size(img_data_raw)

img_FA_fsl_all = nii_FA.img;
img_MD_fsl_all = nii_MD.img;
img_S0_fsl_all = nii_S0.img;

img_V1_fsl_all = nii_V1.img;
img_V2_fsl_all = nii_V2.img;
img_V3_fsl_all = nii_V3.img;

img_L1_fsl_all = nii_L1.img;
img_L2_fsl_all = nii_L2.img;
img_L3_fsl_all = nii_L3.img;

map_JHU_wm_label = JHU_wm_label.img;
map_ATR_R = ATR_R.img;
map_Cigulum = Cigulum.img;
map_CT_R = CT_R.img;
map_Fm_R = Fm_R.img;
map_IFOF_R = IFOF_R.img;
map_SLF_R = SLF_R.img;
map_ROI = ROI.img;

clearvars nii_data nii_FA nii_MD nii_MO nii_S0 nii_V1 nii_V2 nii_V3 nii_L1...
    nii_L2 nii_L3 JHU_wm_label ATR_L Cigulum CT_L Fm_L IFOF_L...
    SLF_L ROI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check crossing fibers in ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_thresh = 10;

sum(sum(sum(map_ATR_L)))
sum(sum(sum(map_Cigulum)))
sum(sum(sum(map_CT_L)))
sum(sum(sum(map_FM_L)))
sum(sum(sum(map_Fm_L)))
sum(sum(sum(map_IFOF_L)))
sum(sum(sum(map_ILF_L)))
sum(sum(sum(map_SLF_L)))
sum(sum(sum(map_ROI)))


sum(sum(sum(map_SLF_L>mask_thresh&map_SCR_L>mask_thresh)))
sum(sum(sum(map_SLF_L>mask_thresh&map_PCR_L>mask_thresh)))
sum(sum(sum(map_SLF_L>mask_thresh&map_BCC_L>mask_thresh)))
sum(sum(sum(map_SCR_L>mask_thresh&map_BCC_L>mask_thresh)))
sum(sum(sum(map_PCR_L>mask_thresh&map_BCC_L>mask_thresh)))

% save(strcat('/Users/hao/Dropbox/DMRI_code/','nii_data_header.mat'));
