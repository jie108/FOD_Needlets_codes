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
% addpath(strcat(path_folder,'mri/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load fitted space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_range = 44:59; %71:105;%71:120; %86:95; %71:105;  %71:110; %120:124;  %100:155    101:160
y_range = 75:90; %41:125;%46:145; %86:90;  %41:130; %137:141;  %87:148     91:180;   
z_range = 55:70; %51:90;%51:100; %71:85; %61:90;  %46:105; %38:40;    %21:47      21:50;
voxel_sub = strcat('x',num2str(x_range(1)),'-',num2str(x_range(end)),'y',num2str(y_range(1)),'-',num2str(y_range(end)),'z',num2str(z_range(1)),'-',num2str(z_range(end)));

b = 5000;
load(strcat(path_folder,'fitting_b',num2str(b),'/','space_index', voxel_sub,'.mat'));

clear options;
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


x_sub = 4:8;
y_sub = 6:14;
z_sub = 3:11;
n1 = length(x_sub);
n2 = length(y_sub);
n3 = length(z_sub);

oax_left = 0;
oax_bottom = 0;
oax_width = 1;
oax_height = 1;

ax_left = 0.05;
ax_bottom = 0.05;
ax_width = 0.7;
ax_height = 0.9;

width_space = ax_width/n2;
height_space = ax_height/n3;  

for i1 = x_sub  %% x-perspective; for z-perspetive, change k1, n1, to k3, n3 and indexing to the third one
        
    figure;
        h    = [];
        h(1) = subplot(2,3,1);
        h(2) = subplot(2,3,2);
        h(3) = subplot(2,3,3);
        h(4) = subplot(2,3,4);
        h(5) = subplot(2,3,5);
        h(6) = subplot(2,3,6);

        %% FA map
        FA_temp_map = dti_data_reindex(data_FA(i1,y_sub,z_sub));
        imagesc(FA_temp_map,'Parent',h(1));
        set(gca,'YDir','normal')
        axis off
        colormap('gray'); %% gray scale 
        caxis manual  %% 
        caxis([FA_bottom FA_top]);  %% fix scale such that it is comparable across slices
%         colorbar;

        %% MD map
        MD_temp_map = dti_data_reindex(data_MD(i1,y_sub,z_sub));
        imagesc(MD_temp_map,'Parent',h(2));
        colormap('gray');
        caxis manual
        caxis([MD_bottom MD_top]);
%         colorbar;

        %% Color leading eigenvector map
        Eig1_FAcorrected_temp = dti_data_reindex(squeeze(Eig1_FAcorrected(i1,y_sub,z_sub,:)));
%         Eig1_fig = figure;
        axis_FA_fig = axes;
        imagesc(Eig1_FAcorrected_temp,'Parent',h(3));
        axis(axis_FA_fig,'off')
        title('FA, MD, Colormap; Top est.; Bottom fsl')

%%%%%%%%%%%%%%%%% fsl results for comparison
        %% FA map
        FA_temp_map_fsl = dti_data_reindex(img_FA_fsl(i1,y_sub,z_sub));
        imagesc(FA_temp_map_fsl,'Parent',h(4));
        colormap('gray');
        caxis manual
        caxis([FA_bottom FA_top]);
%         colorbar;

        %% MD map
        MD_temp_map_fsl = dti_data_reindex(img_MD_fsl(i1,y_sub,z_sub)*1e3);
        imagesc(MD_temp_map_fsl,'Parent',h(5));
        colormap('gray');
        caxis manual
        caxis([MD_bottom MD_top]);
%         colorbar;

        %% Color leading eigenvector map
        Eig1_FAcorrected_temp_fsl = dti_data_reindex(squeeze(Eig1_fsl_FAcorrected(i1,y_sub,z_sub,:)));
%         Eig1_fig = figure;
        axis_FA_fig_fsl = axes;
        imagesc(Eig1_FAcorrected_temp_fsl,'Parent',h(6));
        axis(axis_FA_fig_fsl,'off')

      
        display(i1);
        savefig(strcat(save_path,'sub99_FA_MD_color_x',num2str(x_range(i1)),'.fig'))
        saveas(gcf,strcat(save_path,'figure/sub99_FA_MD_color_x',num2str(x_range(i1)),'.png'));
        close(gcf);
end

for i1 = x_sub
            FA_temp_map = dti_data_reindex(squeeze(data_FA(i1,y_sub,z_sub)));
            %% plot est. FOD on FA background inversely scaled with MD value
            FA_fig = figure('units','normalized','position',[0 0 0.8 1]);
            axis_FA_fig = axes;
            colormap(axis_FA_fig, gray);
            imagesc(FA_temp_map);
%             colorbar;
            caxis(axis_FA_fig,[FA_bottom FA_top]);
            ax = gca;  
            ax.Units = 'normalized';
            ax.OuterPosition = [oax_left, oax_bottom, oax_width, oax_height];
            ax.Position = [ax_left, ax_bottom, ax_width, ax_height];
            axis(axis_FA_fig,'off')

            
            hold on;
            for i2 = y_sub
                for i3 = z_sub
                    if(data_MD(i1,i2,i3)<1.5)
                        fig_factor = 0.2;
                    else
                        fig_factor = 0.2+(data_MD(i1,i2,i3)-1.5)/MD_top;
                    end

                    ax_temp_coordinate = [ax_left+(i2-(y_sub(1)-1)-1)*width_space+width_space*(fig_factor)/2, ax_bottom+(i3-(z_sub(1)-1)-1)*height_space+height_space*(fig_factor)/2, width_space*(1-fig_factor), height_space*(1-fig_factor)];
                    ax_temp = axes('Position', ax_temp_coordinate);
                    axis(ax_temp,'off');
                    colormap(ax_temp,jet);
                    options.use_axis = ax_temp_coordinate;
                    plot_spherical_function(v_p,f_p,squeeze(fod_sele_SN_all(i1,i2,i3,:)),options);
%                     alpha(0.25)
                    lightangle(pi/2,0)
                    view([1 0 0])
                end
            end
            hold off;
            savefig(strcat(save_path,'sub99_SN_x', num2str(x_range(i1)),'.fig'));          
            saveas(gcf,strcat(save_path,'figure/sub99_SN_x',num2str(x_range(i1)),'.png'));
            close(gcf);
            display(i1);
%fod_sele_SN_all, fod_SH_all8, fod_SCSD_lmax8_all, fod_SCSD_lmax12_all,
%SN_x, SH8_x, SCSD8_x, SCSD12_x,
end


for i1 = x_sub  %% x-perspective; for z-perspetive, change k1, n1, to k3, n3 and indexing to the third one
        
        %% FA map
        figure;
        FA_temp_map = dti_data_reindex(data_FA(i1,y_sub,z_sub));
        imagesc(FA_temp_map);
%         set(gca,'YDir','normal')
%         axis off
        colormap('gray'); %% gray scale 
        caxis manual  %% 
        caxis([FA_bottom FA_top]);  %% fix scale such that it is comparable across slices
        colorbar;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 11 10]); %x_width=11cm y_width=10cm
        saveas(gcf,strcat(save_path,'figure/sub99_FA_x',num2str(x_range(i1)),'.png'));
end
%{
for i1 = x_sub  %% x-perspective; for z-perspetive, change k1, n1, to k3, n3 and indexing to the third one
         %% MD map
    figure;
         %% MD map
        MD_temp_map = dti_data_reindex(data_MD(i1,y_sub,z_sub));
        imagesc(MD_temp_map);
        colormap('gray');
        caxis manual
        caxis([MD_bottom MD_top]);
%         colorbar;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 11 10]); %x_width=11cm y_width=10cm
%         saveas(gcf,strcat(save_path,'figure/MD_x',num2str(x_range(i1)),'.png'));
end
%}
for i1 = x_sub  %% x-perspective; for z-perspetive, change k1, n1, to k3, n3 and indexing to the third one
        
    figure;
          %% Color leading eigenvector map
        Eig1_FAcorrected_temp = dti_data_reindex(squeeze(Eig1_FAcorrected(i1,y_sub,z_sub,:)));
%         Eig1_fig = figure;
        axis_FA_fig = axes;
        imagesc(Eig1_FAcorrected_temp);
        axis(axis_FA_fig,'off')
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 10 10]); %x_width=11cm y_width=10cm
        saveas(gcf,strcat(save_path,'figure/sub99_color_x',num2str(x_range(i1)),'.png'));
end

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% new figures FA, MD map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i1 = x_sub
    fig_temp = figure;
    fig_temp.Units = 'normalized';
    set(fig_temp, 'Position', [0 0 1 1]); 
    h1 = subplot(1,4,1);
    %% Brain region
    brain_temp = imread(strcat(save_path,'/figure/','fsl_x',num2str(x_range(i1)),'.png'));
    imshow(brain_temp);
%     title('Region')
    
    h2 = subplot(1,4,2);
    %% Color leading eigenvector map
    Eig1_FAcorrected_temp = dti_data_reindex(squeeze(Eig1_FAcorrected(i1,y_sub,z_sub,:)));
    axis_FA_fig = axes;
    imagesc(Eig1_FAcorrected_temp,'Parent',h2);
    axis(axis_FA_fig,'off')
%     title(h2,'Colormap')
    
    h3 = subplot(1,4,3);
    originalSize3 = get(gca, 'Position');

    h4 = subplot(1,4,4);
    originalSize4 = get(gca, 'Position');

    %% FA map
    FA_temp_map = dti_data_reindex(data_FA(i3,y_sub,z_sub));
    imagesc(FA_temp_map,'Parent',h3);
%     title(h3,'FA');
    colormap('gray');
    set(h3,'yticklabel',[])
    cb3 = colorbar(h3);
    caxis(h3,[FA_bottom FA_top]);
    
    
    %% MD map
    MD_temp_map = dti_data_reindex(data_MD(i3,y_sub,z_sub));
    imagesc(MD_temp_map,'Parent',h4);
    colormap('gray');
    set(h4,'xticklabel',[])
%     title(h4,'MD');
    colormap('gray');
    set(h4,'yticklabel',[])
    cb4 = colorbar(h4);
    caxis(h4,[MD_bottom MD_top]);
    
%     set(h3, 'Position', originalSize3); 
%     set(h4, 'Position', originalSize4); 
    saveas(gcf,strcat(save_path,'figure/','fsl_color_FA_MD_x',num2str(x_range(i1)),'.png'));
    close(gcf);
end
%}