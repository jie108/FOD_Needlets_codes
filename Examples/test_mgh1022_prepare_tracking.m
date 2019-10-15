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
path_folder = strcat(path_save, subject_id, '/fitting/');

addpath(path_folder);

load(strcat(path_folder,'space_indexx108-123y124-139z37-42/','space_indexx108-123y124-139z37-42.mat'));
% load(strcat(savepath,'for_tracking.mat'));

savepath = strcat(path_folder,'space_indexx108-123y124-139z37-42/');

%%%vertex construction options 
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;

%%%plotting options 
options.spherical = 1;
% options for the display
options.use_color = 1;
options.color = 'wavelets';
options.use_elevation = 2;
options.rho = 0.5;
options.scaling = 1.5;
% for draw_fiber
plot_rho = options.rho;
plot_sacle = options.scaling;

% denser grid for interpolation and plots 
[v_p,f_p] = compute_semiregular_sphere(5,options);
pos_p = v_p{end};
phi_p = atan2(pos_p(2,:),pos_p(1,:))/(2*pi);   %%phi: azimuthal  angle, [0,2)
phi_p = phi_p+(phi_p<0);
theta_p = acos(pos_p(3,:))/(pi);             %% theta: polar angle, [0,1)


xgrid_sp = 1.3672; 
ygrid_sp = 1.3672;
zgrid_sp = 2.7;

braingrid = zeros(3,n1,n2,n3);
for i = 1:n1
  for j = 1:n2
    for k = 1:n3
      braingrid(:,i,j,k) = [(i-n1/2-.5)*xgrid_sp (j-n2/2-.5)*ygrid_sp (k-n3/2-.5)*zgrid_sp];
    end
  end
end

kmin = 40;
peak_thresh = 0.25;
num_fib_cut = 3;
Dis = squareform(pdist(pos_plot','cosine'));


n_fiber = zeros(1,n1*n2*n3);
rmap = zeros(1,n1*n2*n3);

vec = [];
loc = [];
map = [];

for k=1:n3
    for j=1:n2
        for i=1:n1
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            [~, ~, ~, ~, ~, peak_pos_SN_final_RSSdiff] = FOD_peak(fod_sele_SN_RSSdiff(i,j,k,:), Dis, kmin, peak_thresh, pos_p, theta_p, phi_p);
           
            n_fiber((k-1)*n1*n2+(j-1)*n1+i) = size(peak_pos_SN_final_RSSdiff,2);
            if((size(peak_pos_SN_final_RSSdiff,2)>=1)&&(size(peak_pos_SN_final_RSSdiff,2)<=num_fib_cut)) %&&(size(peak_pos_SN_final_RSSdiff,2)<=4))
                vec = [vec, peak_pos_SN_final_RSSdiff];
                map = [map repmat((k-1)*n1*n2+(j-1)*n1+i,1,size(peak_pos_SN_final_RSSdiff,2))];
                loc = [loc repmat(braingrid(:,i,j,k),1,size(peak_pos_SN_final_RSSdiff,2))];
            else
                vec = [vec, [NaN; NaN; NaN]];
                map = [map repmat((k-1)*n1*n2+(j-1)*n1+i,1,1)];
                loc = [loc repmat(braingrid(:,i,j,k),1,1)];                
            end
            %map = [map repmat((k-1)*n1*n2+(j-1)*n1+i,1,size(peak_pos_SN_final_RSSdiff,2))];
        end
        display(j);
    end
    display(k);
end

n_fiber2 = [];
n_fiber(n_fiber>num_fib_cut)=0;
for i=1:n1*n2*n3
    temp_v = find(map==i);
    rmap(i) = temp_v(1);
    n_fiber2 = [n_fiber2 repmat(n_fiber(i),1,max(n_fiber(i),1))];
end

vec = vec';
loc = loc';
save(strcat(savepath,'for_tracking_cut',num2str(num_fib_cut),'.mat'),'vec','loc','n_fiber','n_fiber2','map'...
    ,'rmap','braingrid','xgrid_sp','ygrid_sp','zgrid_sp','n1','n2','n3','peak_thresh','num_fib_cut');
