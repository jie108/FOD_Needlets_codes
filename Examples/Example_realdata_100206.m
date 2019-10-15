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
% subject_id = '100206';
% path_folder = strcat(path_save, subject_id, '/T1w/', 'Diffusion/');
subject_id = 'mgh_1011';
path_folder = strcat(path_save, subject_id, '/diff/preproc/');

addpath(path_folder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load .nii data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load nii data (raw)
bvec = load('bvecs_moco_norm.txt'); 
bval = load('bvals.txt');
% nii_data = load_untouch_nii(strcat(path_folder,'mri/','diff_preproc_new.nii.gz'));
% nii_data_old = load_nii(strcat(path_folder,'mri/','diff_preproc.nii.gz'));

% reslice_nii(strcat(path_folder,'mri/','diff_preproc.nii.gz'),strcat(path_folder,'mri/','diff_preproc_new.nii.gz'));
% reslice_nii(strcat(path_folder,'mri/','dti_FA.nii.gz'),strcat(path_folder,'mri/','dti_FA_new.nii.gz'));
% reslice_nii(strcat(path_folder,'mri/','dti_MO.nii.gz'),strcat(path_folder,'mri/','dti_MO_new.nii.gz'));
% reslice_nii(strcat(path_folder,'mri/','dti_MD.nii.gz'),strcat(path_folder,'mri/','dti_MD_new.nii.gz'));
% reslice_nii(strcat(path_folder,'mri/','dti_S0.nii.gz'),strcat(path_folder,'mri/','dti_S0_new.nii.gz'));
% reslice_nii(strcat(path_folder,'mri/','dti_V1.nii.gz'),strcat(path_folder,'mri/','dti_V1_new.nii.gz'));
% reslice_nii(strcat(path_folder,'mri/','dti_V2.nii.gz'),strcat(path_folder,'mri/','dti_V2_new.nii.gz'));
% reslice_nii(strcat(path_folder,'mri/','dti_V3.nii.gz'),strcat(path_folder,'mri/','dti_V3_new.nii.gz'));
% reslice_nii(strcat(path_folder,'mri/','dti_L1.nii.gz'),strcat(path_folder,'mri/','dti_L1_new.nii.gz'));
% reslice_nii(strcat(path_folder,'mri/','dti_L2.nii.gz'),strcat(path_folder,'mri/','dti_L2_new.nii.gz'));
% reslice_nii(strcat(path_folder,'mri/','dti_L3.nii.gz'),strcat(path_folder,'mri/','dti_L3_new.nii.gz'));

nii_data = load_nii(strcat(path_folder,'mri/','diff_preproc_new.nii.gz'));

%% load fsl processed data
nii_FA = load_nii('dti_FA_new.nii.gz');
nii_MO = load_nii('dti_MO_new.nii.gz');
nii_MD = load_nii('dti_MD_new.nii.gz');
nii_S0 = load_nii('dti_S0_new.nii.gz');

nii_V1 = load_nii('dti_V1_new.nii.gz');
nii_V2 = load_nii('dti_V2_new.nii.gz');
nii_V3 = load_nii('dti_V3_new.nii.gz');

nii_L1 = load_nii('dti_L1_new.nii.gz');
nii_L2 = load_nii('dti_L2_new.nii.gz');
nii_L3 = load_nii('dti_L3_new.nii.gz');

%% transfer to matrix
img_data_raw = nii_data.img;
size(img_data)

img_FA_fsl_all = nii_FA.img;
img_MD_fsl_all = nii_MD.img;
img_S0_fsl_all = nii_S0.img;

img_V1_fsl_all = nii_V1.img;
img_V2_fsl_all = nii_V2.img;
img_V3_fsl_all = nii_V3.img;

img_L1_fsl_all = nii_L1.img;
img_L2_fsl_all = nii_L2.img;
img_L3_fsl_all = nii_L3.img;

b = 1000;
b_sample = find(bval==b);
img_data_all = img_data_raw(:,:,:,b_sample);
b0_sample = find(bval==0);
img_b0_all = img_data_raw(:,:,:,b0_sample);

clearvars nii_data nii_FA nii_MD nii_MO nii_S0 nii_V1 nii_V2 nii_V3 nii_L1 nii_L2 nii_L3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% real data range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_range = 71:105;  %71:110; %120:124;  %100:155    101:160
y_range = 71:90;  %41:130; %137:141;  %87:148     91:180;   
z_range = 61:90;  %46:105; %38:40;    %21:47      21:50;

img_data = img_data_all(x_range,y_range,z_range,:);
img_b0 = img_b0_all(x_range,y_range,z_range,:);

img_S0_fsl = img_S0_fsl_all(x_range,y_range,z_range);
img_MD_fsl = img_MD_fsl_all(x_range,y_range,z_range);
img_FA_fsl = img_FA_fsl_all(x_range,y_range,z_range);
img_V1_fsl = img_V1_fsl_all(x_range,y_range,z_range,:);
img_L1_fsl = img_L1_fsl_all(x_range,y_range,z_range,:);

[n1, n2, n3, n4] = size(img_data);


%% MLE estimation of S0 and sigma in each voxel

Sigma_mle = zeros(n1,n2,n3);
S0_mle = zeros(n1,n2,n3);
% Est_sigma_all = [];
% Est_S0_all = [];
for k1=1:n1
    for k2=1:n2
        for k3=1:n3
            
            y = squeeze(img_b0(k1,k2,k3,:));
            if(min(y)>0&&max(abs(y))~=Inf)
                options = optimoptions('fminunc','Algorithm','quasi-newton'); %,'TolFun',1e-20);
                options.Display = 'off';
                x0 = double([var(y),mean(y)]);
%                 x00 = double([3600 4000]);
                f = @(x) parameterfun(x,y);
                [x, fval, exitflag,grad] = fminunc(f,x0,options);
%                 [x, fval, exitflag, output] = fminunc(f,x0,options);
%                 Est_sigma_all(i,j,k) = sqrt(x(1));
%                 Est_S0_all(i,j,k) = x(2);
                Sigma_mle(k1,k2,k3) = sqrt(x(1));
                S0_mle(k1,k2,k3) = (x(2));
            end
        end
    end
   	display(k1);
end

Sigma_var = zeros(n1,n2,n3);
S0_mean = zeros(n1,n2,n3);
for k1=1:n1
    for k2=1:n2
        for k3=1:n3
            y = squeeze(img_b0(k1,k2,k3,:));
            Sigma_var(k1,k2,k3) = sqrt(var(y));
            S0_mean(k1,k2,k3) = mean(y);
        end
    end
end


figure
subplot(4,2,1)  % histogram of MLE est. sigma
hist(reshape(Est_sigma_all,1,k1*k2*k3))
title('Sigma')
subplot(4,2,2)  % histogram of MLE est. S0
hist(reshape(Est_S0_all,1,k1*k2*k3))
title('S0')
subplot(4,2,3)  % MLE sigma vs. Variance est. sigma
scatter(reshape(Est_sigma_all,1,k1*k2*k3),reshape(Est_sigma_var_all,1,k1*k2*k3))
title('var vs. Sigma')
subplot(4,2,4)  % MLE S0 vs. Mean est. S0
scatter(reshape(Est_S0_all,1,k1*k2*k3),reshape(Est_S0_mean_all,1,k1*k2*k3))
title('mean vs S0')
subplot(4,2,5)  % MLE S0 vs. MLE sigma
scatter(reshape(Est_S0_all,1,k1*k2*k3), reshape(Est_sigma_all,1,k1*k2*k3))
title('Sigma vs. S0')
subplot(4,2,6)  % fsl MD vs MLE S0
scatter( reshape(img_MD_fsl_temp,1,k1*k2*k3),reshape(Est_S0_all,1,k1*k2*k3))
title('S0 vs. MD')
subplot(4,2,7)  % fsl MD vs MLE sigma
scatter( reshape(img_MD_fsl_temp,1,k1*k2*k3),reshape(Est_sigma_all,1,k1*k2*k3))
title('Sigma vs. MD')
subplot(4,2,8) % SNR as ratio of MLE S0 over MLE sigma
hist(reshape(Est_S0_all./Est_sigma_all,1,k1*k2*k3),20)
title('SNR')

%% SNR statistics
% median(reshape(Est_S0_all,1,k1*k2*k3))/median(reshape(Est_sigma_all,1,k1*k2*k3))
Est_S0_all_vec = reshape(Est_S0_all,1,k1*k2*k3);
Est_sigma_all_vec = reshape(Est_sigma_all,1,k1*k2*k3);

SNR = Est_S0_all./Est_sigma_all;
SNR_vec = reshape(SNR,1,k1*k2*k3);

% SNR_vec((isnan(SNR_vec))) = 0; 

min(SNR_vec(SNR_vec>0))
max(SNR_vec(SNR_vec>0))
mean(SNR_vec(SNR_vec>0))
median(SNR_vec(SNR_vec>0))

quantile(SNR_vec,0.05)
quantile(SNR_vec,0.95)
figure
hist(SNR_vec,50)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Single tensor model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = k1*k2*k3;
n = 41;
S_min = min(min(img_temp(img_temp>0)));
S_max = max(max(img_temp(img_temp<Inf)));

S_all_temp = zeros(41,N);
for i = 1:41
%     S_all_temp(i,:) = reshape(max (S_min,img_temp(:,:,:,i+5)),1,N);
%     S_all_temp(i,:) = reshape(min(S_max,img_temp(:,:,:,i+5)),1,N);
    S_all_temp(i,:) = reshape(img_temp(:,:,:,i+5),1,N);
    S_all_temp(i,:) = S_all_temp(i,:)./reshape(Est_S0_all,1,N);  %% use S0_fsl
end

index_negative = [];
index_infinity = [];
for i=1:N
    if(min(S_all_temp(:,i))<0)
        index_negative = [index_negative i];
    end
    if(sum(abs(S_all_temp(:,i)==Inf))>0)
        index_infinity = [index_infinity i];
    end
end
index_temp = setdiff(1:N,union(index_negative,index_infinity));




X = zeros(n,6);
X(:,1) = bvec(1,6:46).^2;
X(:,2) = bvec(2,6:46).^2;
X(:,3) = bvec(3,6:46).^2;
X(:,4) = 2*bvec(1,6:46).*bvec(2,6:46);
X(:,5) = 2*bvec(1,6:46).*bvec(3,6:46);
X(:,6) = 2*bvec(2,6:46).*bvec(3,6:46);

% N = length(index_temp);

FA_temp = zeros(1,N);
eval_temp = zeros(3,N);
evec_temp = zeros(3,3,N);

%% single tensor model
for i=1:N   %length(index_temp)
    %i = index_temp(j);
    l_S = log(S_all_temp(:,i));
    D_est_temp = -inv((X'*X))*X'*l_S;
    [D_nl, iter, DWI_est] = LM_dti(X,S_all_temp(:,i),D_est_temp,1e-10);
    D_matrix_nl = [D_nl(1) D_nl(4) D_nl(5);D_nl(4) D_nl(2) D_nl(6);D_nl(5) D_nl(6) D_nl(3)];
    if(sum(isnan(D_matrix_nl))==0)
      [egvec_nl, egval_nl] = eig(D_matrix_nl);
    end
    eval_temp(:,i) = diag(egval_nl);
    evec_temp(:,:,i) = egvec_nl;
    FA_temp(i) = sqrt(1/2)*sqrt(((eval_temp(1,i)-eval_temp(2,i))^2+(eval_temp(1,i)-eval_temp(3,i))^2+(eval_temp(3,i)-eval_temp(2,i))^2))/sqrt((eval_temp(1,i)^2+eval_temp(2,i)^2+eval_temp(3,i)^2)); 
end


MD_temp = sum(eval_temp,1)./3;
eval_ratio23_temp = eval_temp(2,:)./eval_temp(1,:);
eval_ratio_temp = eval_temp(3,:)*2./(eval_temp(1,:)+eval_temp(2,:));
index_ttemp = find(FA_temp<1&MD_temp>=0);
data_FA = reshape(FA_temp,k1,k2,k3);
data_MD = reshape(MD_temp,k1,k2,k3);

FA_show = FA_temp(index_ttemp);
MD_show = MD_temp(index_ttemp);
eval_ratio_show = eval_ratio_temp(index_ttemp);
eval3_show = eval_temp(3,index_ttemp);
eval_show = eval_temp(:,index_ttemp);
eval_ratio23_show = eval_ratio23_temp(index_ttemp);
SNR_show = reshape(Est_S0_all./Est_sigma_all,1,k1*k2*k3);
median(SNR_show(~isnan(SNR_show)))
mean(SNR_show(~isnan(SNR_show)))

idx_response = find(FA_show<1&FA_show>0.8&eval_show(1,:)>0&eval_show(2,:)>0&eval_show(3,:)>0&eval_ratio23_show<1.5);
size(idx_response)

temp = sort(eval_ratio_show(idx_response),'descend');
ratio_response = temp(floor(size(temp,2)/2));

id_temp = find(eval_ratio_show==ratio_response);
b_factor = median(eval_show(3,id_temp));


figure % compare single tensor model estimation vs fsl estimation
subplot(2,1,1)
scatter(FA_show, reshape(img_FA_fsl_temp,1,numel(img_FA_fsl_temp)))
subplot(2,1,2)
scatter(MD_show, reshape(img_MD_fsl_temp,1,numel(img_MD_fsl_temp)))


figure 
subplot(2,3,1)
hist(reshape(Est_S0_all,1,numel(Est_S0_all)),20)
set(gca, 'FontSize', 15)
title('S0','FontSize',25)
subplot(2,3,2)
hist(reshape(log(Est_sigma_all),1,numel(Est_sigma_all)),20)
set(gca, 'FontSize', 15)
title('log Sigma','FontSize',25)
subplot(2,3,3)
hist(reshape(Est_S0_all./Est_sigma_all,1,k1*k2*k3),20)
set(gca, 'FontSize', 15)
title('SNR','FontSize',25)
subplot(2,3,4)
hist(FA_show,25)
set(gca, 'FontSize', 15)
title('FA','FontSize',25)
subplot(2,3,5)
hist(MD_show,20)
set(gca, 'FontSize', 15)
title('MD','FontSize',25)
subplot(2,3,6)
scatter(FA_show,MD_show)
set(gca, 'FontSize', 15)
title('MD vs. FA','FontSize',25)
% subplot(2,4,6)
% hist(eval3_show,20)
% title('Eval. 3','FontSize',15)
% subplot(2,4,7)
% Ctrs = [0:2:20];
% Xtrs = hist(eval_ratio_show,Ctrs);
% bar(Ctrs, Xtrs)
% title('Eval. Ratio','FontSize',15)

% eval_ratio23_fsl = img_L2_temp./img_L3_temp;
% eval_ratio_fsl = 2*img_L1_temp./(img_L3_temp+img_L2_temp);

%% set uniform color range for heat map later
FA_temp_restrict = min(1,FA_temp);
MD_temp_restrict = max(0,MD_temp);
FA_top = max(FA_temp_restrict);
FA_bottom = min(FA_temp_restrict);
MD_top = max(MD_temp_restrict);
MD_bottom = min(MD_temp_restrict);

% for colormap
Eig1_FAcorrected = zeros(k1,k2,k3,3);
Eig1_fsl_FAcorrected = zeros(k1,k2,k3,3);

for i = 1:N
    [i1,i2,i3] = ind2sub([k1,k2,k3],i);
    Eig1_FAcorrected(i1,i2,i3,:) = abs(evec_temp(:,3,i)).*FA_temp(i);
    Eig1_fsl_FAcorrected(i1,i2,i3,:) = abs(img_V1_fsl_temp(i1,i2,i3,:)).*img_FA_fsl_temp(i1,i2,i3);
end

%%%%%%%%%%%%%%%%%%% FA, MD and color maps
%{
    %% 
    for k = 1:k3  %% x-perspective; for z-perspetive, change k1, n1, to k3, n3 and indexing to the third one
        figure;
        h    = [];
        h(1) = subplot(2,3,1);
        h(2) = subplot(2,3,2);
        h(3) = subplot(2,3,3);
        h(4) = subplot(2,3,4);
        h(5) = subplot(2,3,5);
        h(6) = subplot(2,3,6);

        %% FA map
        FA_temp_map = dti_data_reindex(data_FA(:,:,k));
        imagesc(FA_temp_map,'Parent',h(1));
        set(gca,'YDir','normal')
        axis off
        colormap('gray'); %% gray scale 
        caxis manual  %% 
        caxis([FA_bottom FA_top]);  %% fix scale such that it is comparable across slices
%         colorbar;

        %% MD map
        MD_temp_map = dti_data_reindex(data_MD(:,:,k));
        imagesc(MD_temp_map,'Parent',h(2));
        colormap('gray');
        caxis manual
        caxis([MD_bottom MD_top]);
%         colorbar;

        %% Color leading eigenvector map
        Eig1_FAcorrected_temp = dti_data_reindex(squeeze(Eig1_FAcorrected(:,:,k,:)));
%         Eig1_fig = figure;
        axis_FA_fig = axes;
        imagesc(Eig1_FAcorrected_temp,'Parent',h(3));
        axis(axis_FA_fig,'off')
        title('FA, MD, Colormap; Top est.; Bottom fsl')

%%%%%%%%%%%%%%%%% fsl results for comparison
        %% FA map
        FA_temp_map_fsl = dti_data_reindex(img_FA_fsl_temp(:,:,k));
        imagesc(FA_temp_map_fsl,'Parent',h(4));
        colormap('gray');
        caxis manual
        caxis([FA_bottom FA_top]);
%         colorbar;

        %% MD map
        MD_temp_map_fsl = dti_data_reindex(img_MD_fsl_temp(:,:,k));
        imagesc(MD_temp_map_fsl,'Parent',h(5));
        colormap('gray');
        caxis manual
        caxis([MD_bottom MD_top]);
%         colorbar;

        %% Color leading eigenvector map
        Eig1_FAcorrected_temp_fsl = dti_data_reindex(squeeze(Eig1_fsl_FAcorrected(:,:,k,:)));
%         Eig1_fig = figure;
        axis_FA_fig_fsl = axes;
        imagesc(Eig1_FAcorrected_temp_fsl,'Parent',h(6));
        axis(axis_FA_fig_fsl,'off')

        display(k);
    end


    options.use_axis = 0;
%}




clear options


