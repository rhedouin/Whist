% Create figure from brain for WM paper
clear
close all

base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/';
orients9_fa20_folder = [base_folder 'concatenate_signals_9_orientations/parameter_maps/noise4/fa-20/'];
orients3_fa20_mean_folder = [base_folder 'concatenate_signals_3_orientations/mean_fa-20/'];

parameter_list = {'FVF', 'gRatio', 'R2Myelin', 'R2IntraExtraAxonal', 'weight', 'xiMyelin'};

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_mask_coronal_73.nii.gz'));

fa = 'fa-20';

for k = 1:length(parameter_list)
    
    parameter = parameter_list{k};
    parameter_map_orients9.(parameter) = load_nii_img_only([orients9_fa20_folder 'BrainSample-2_ses-03_' parameter '_fa-20_20_directions_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    parameter_map_orients9.(parameter)(mask==0) = 10000;

    parameter_map_orients9_permute.(parameter) = zeros(82, 128, 128);
    parameter_map_orients9_permute.(parameter)(1:10, :, :) = 10000;
    parameter_map_orients9_permute.(parameter)(11:end, :, :) = permute(parameter_map_orients9.(parameter), [3 2 1]);
    parameter_map_orients9_permute.(parameter) = parameter_map_orients9_permute.(parameter)(end:-1:1, end:-1:1, :);

    
    parameter_map_orients3.(parameter) = load_nii_img_only([orients3_fa20_mean_folder 'BrainSample-2_ses-03_' parameter '_mean_fa-20_20_directions_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    parameter_map_orients3.(parameter)(mask==0) = 0;

    parameter_map_orients3_permute.(parameter) = zeros(82, 128, 128);
    parameter_map_orients3_permute.(parameter)(1:10, :, :) = 0;
    parameter_map_orients3_permute.(parameter)(11:end, :, :) = permute(parameter_map_orients3.(parameter), [3 2 1]);
    parameter_map_orients3_permute.(parameter) = parameter_map_orients3_permute.(parameter)(end:-1:1, end:-1:1, :);

    diff.(parameter) = abs(parameter_map_orients9_permute.(parameter) - parameter_map_orients3_permute.(parameter));
end

xslice = 0;
yslice = 128-72;
zslice = 0;

edge_x = 0;
edge_y = 0;
edge_z = 10;

for k = 1:length(parameter_list)
    parameter = parameter_list{k};
    
%     figure('Name',parameter)
%     imagesc(squeeze(parameter_map_orients9_permute.(parameter)(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));
%     colormap('gray');
%     cb(k) = colorbar('south');

% %     cb(k) = colorbar('southoutside');
%  
%     x1=get(gca,'position');
%     x=get(cb(k),'Position');
% 
%     x(1) = 0.275;
%     x(2) = 0.075;
%     x(3) = 0.5;
%     x(4) = 0.03;
%      
%     set(cb(k),'Position',x)
%     set(gca,'position',x1)
% 
%     if k == 1
%         caxis([0 0.8])
%     elseif k == 2
%         caxis([0.5 0.85])
%     elseif k == 3
%         caxis([30 120])
%     elseif k == 4
%         caxis([10 60])
%     elseif k == 5
%         caxis([0 3])
%     elseif k == 6
%         caxis([-0.2 0.2])
%     end
    
%     figure('Name',parameter)
%     imagesc(squeeze(parameter_map_orients3_permute.(parameter)(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));
%     colormap('gray');
%     cb(k) = colorbar('south');
%     caxis([30 120])
        
    figure('Name',parameter)
    imagesc(squeeze(diff.(parameter)(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));
    colormap('gray');
    cb(k) = colorbar('south');
   
    if k == 1
        caxis([0 0.4])
    elseif k == 2
        caxis([0 0.175])
    elseif k == 3
        caxis([0 40])
    elseif k == 4
        caxis([0 25])
    elseif k == 5
        caxis([0 1])
    elseif k == 6
        caxis([0 0.2])
    end
    
        colormap('gray');
    cb(k) = colorbar('south');
 
    x1=get(gca,'position');
    x=get(cb(k),'Position');

    x(1) = 0.275;
    x(2) = 0.075;
    x(3) = 0.5;
    x(4) = 0.03;
     
    set(cb(k),'Position',x)
    set(gca,'position',x1)

    axis off
    set(gca, 'FontSize', 20)
    
end
return;








