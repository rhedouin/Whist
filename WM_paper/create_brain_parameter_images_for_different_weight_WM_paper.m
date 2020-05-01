% Create figure from brain for WM paper
clear
close all

concatenate_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals_9_orientations/';
parameter_folder = [concatenate_folder 'parameter_maps/noise4/'];

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_gre_orientation-4_brain_mask_all_register_coronal.nii.gz'));

fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60'}

for k = 1:length(fa_list)
    fa = fa_list{k};

    flip_angle_folder = [parameter_folder fa '/'];

    toto = load_untouch_nii([flip_angle_folder 'BrainSample-2_ses-03_weight_' fa '_20_directions_polyfit_cartesian_with_theta_noise4.nii.gz']);
    weight_parameter_map{k} = toto.img;
%     weight_parameter_map{k} = load_untouch_nii([flip_angle_folder 'BrainSample-2_ses-03_weight_' fa '_20_directions_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
%     weight_parameter_map{k} = load_nii_img_only([flip_angle_folder 'BrainSample-2_ses-03_weight_' fa '_20_directions_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
%     weight_parameter_map{k}(mask==0) = 10000;
    
    weight_parameter_map_permute{k} = permute(weight_parameter_map{k}, [3 2 1]);
    weight_parameter_map_permute{k} = weight_parameter_map_permute{k}(end:-1:1, end:-1:1, :);

    toto.img =   weight_parameter_map_permute{k};
    save_untouch_nii(toto, [flip_angle_folder 'BrainSample-2_ses-03_weight_' fa '_20_directions_polyfit_cartesian_with_theta_noise4_registered_test.nii.gz']);

    return;
end



xslice = 0;
yslice = 128 - 65;
zslice = 0;

edge_x = 0;
edge_y = 0;
edge_z = 20;



for k = 1:length(fa_list) 
    fa = fa_list{k};

    figure('Name',fa)
    imagesc(squeeze(weight_parameter_map_permute{k}(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));

    colormap('gray');
    cb(k) = colorbar('south');

    x=get(cb(k),'Position');
    x(1) = 0.2;
    x(2) = 0.01;
    x(3) = 0.6;
    x(4) = 0.03;

    set(cb(k),'Position',x)

    caxis([0.5 2.5])

    
    axis off
    set(gca, 'FontSize', 20)
end




