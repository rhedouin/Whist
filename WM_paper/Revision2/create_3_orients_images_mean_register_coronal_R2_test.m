clear;
close all


base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/';
concatenate_folder = [base_folder 'concatenate_signals_3_orientations/'];
experiment_name = 'BrainSample2LorentzianCorrection';

fa = 'fa-35';
parameter = 'FVF';
list_orientations = {'457', '157', '246', '245', '679', '269', '467', '579', '125', '234'};

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_mp2rage_orientation-4_t1_brain_2_BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1_brain_register_masked_mask_f05.nii.gz'));


xslice = 0;
yslice = 128-72;
zslice = 0;

edge_x = 0;
edge_y = 0;
edge_z = 10;

parameter_map_permute_total = zeros(82, 128, 128);

%%%%%%%%%%%%%% Create R2 maps
for m = 1 : length(list_orientations)
    current_orientation = list_orientations{m};
    orientation_folder = [concatenate_folder current_orientation 'orientations/parameter_maps/' experiment_name '/'];
    
    fa_folder = [orientation_folder fa '/'];
    
    input_filename = [fa_folder 'BrainSample-2_ses-03_' parameter '_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
    
    parameter_map{m} = load_nii_img_only([fa_folder 'BrainSample-2_ses-03_' parameter '_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);

    parameter_map{m}(mask==0) = 10000;

    parameter_map_permute{m} = zeros(82, 128, 128);
    parameter_map_permute{m}(1:10, :, :) = 10000;
    parameter_map_permute{m}(11:end, :, :) = permute(parameter_map{m}, [3 2 1]);
    parameter_map_permute{m} = parameter_map_permute{m}(end:-1:1, end:-1:1, :);
    
    parameter_map_permute_total = parameter_map_permute_total + parameter_map_permute{m};
    
    figure
    imagesc(squeeze(parameter_map_permute{m}(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));
    colormap('gray');
    cb(m) = colorbar('south');
        
    x1=get(gca,'position');
    x=get(cb(m),'Position');
    
    x(1) = 0.275;
    x(2) = 0.075;
    x(3) = 0.5;
    x(4) = 0.03;
    
    set(cb(m),'Position',x)
    set(gca,'position',x1)
    
    caxis([0 0.8])
end

parameter_map_permute_total = parameter_map_permute_total / 10;


figure
imagesc(squeeze(parameter_map_permute{m}(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));
colormap('gray');
cb(m) = colorbar('south');

x1=get(gca,'position');
x=get(cb(m),'Position');

x(1) = 0.275;
x(2) = 0.075;
x(3) = 0.5;
x(4) = 0.03;

set(cb(m),'Position',x)
set(gca,'position',x1)

caxis([0 0.8])



