% Create figure from brain for WM paper
clear
close all

base_folder = ['/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/'];
ref_folder = [base_folder 'references/'];
concatenate_folder = [base_folder 'concatenate_signals_9_orientations/'];
parameter_folder = [concatenate_folder 'parameter_maps/noise4/'];

parameter_list = {'FVF', 'gRatio', 'R2Myelin', 'R2IntraExtraAxonal', 'weight', 'xiMyelin'};

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_mask_coronal_73.nii.gz'));

fa = 'std';
flip_angle_folder = [parameter_folder fa '/'];

for k = 1:length(parameter_list)
    parameter = parameter_list{k};
    parameter_map.(parameter) = load_nii_img_only([flip_angle_folder 'BrainSample-2_ses-03_' parameter '_' fa '_20_directions_polyfit_cartesian_with_theta_noise4_register.nii.gz']);

    parameter_map.(parameter)(mask==0) = 10000;

    parameter_map_permute.(parameter) = zeros(82, 128, 128);
    parameter_map_permute.(parameter)(1:10, :, :) = 10000;
    parameter_map_permute.(parameter)(11:end, :, :) = permute(parameter_map.(parameter), [3 2 1]);
    parameter_map_permute.(parameter) = parameter_map_permute.(parameter)(end:-1:1, end:-1:1, :);
end


xslice = 0;
yslice = 128-72;
zslice = 0;

edge_x = 0;
edge_y = 0;
edge_z = 10;

for k = 1:length(parameter_list)
    parameter = parameter_list{k};
    
    figure('Name',parameter)

    imagesc(squeeze(parameter_map_permute.(parameter)(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));
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

    if k == 1
        caxis([0 0.4])
    elseif k == 2
        caxis([0 0.175])
    elseif k == 3
        caxis([0 60])
    elseif k == 4
        caxis([0 25])
    elseif k == 5
        caxis([0 1.5])
    elseif k == 6
        caxis([0 0.2])
    end
    
    axis off
    %     title(parameter)
    set(gca, 'FontSize', 20)

end











