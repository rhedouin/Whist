% Compute mean and std parameter map for WM paper

parameter_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals_9_orientations/parameter_maps/noise4/';
mean_folder = [parameter_folder 'mean/'];
std_folder = [parameter_folder 'std/'];

fa_list = {'05', '10', '15', '20', '35', '60'}

% parameter_list = {'FVF', 'gRatio', 'T2Myelin', 'T2IntraExtraAxonal', 'weight', 'xiMyelin'};

parameter_map_mean = 0;
parameter_map_example = load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals_9_orientations/parameter_maps/noise4/fa-05/BrainSample-2_ses-03_FVF_fa-05_20_directions_polyfit_cartesian_with_theta_noise4.nii.gz');
dims = size(parameter_map_example);

mask_map = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_gre_orientation-4_brain_mask_all_register.nii.gz'));
mask_erode = mask_map;

 parameter_list = {'T2Myelin', 'T2IntraExtraAxonal'};
output_parameter_list = {'R2Myelin', 'R2IntraExtraAxonal'};

for k = 1:length(parameter_list)
    parameter_map_all = zeros([dims(1) dims(2) dims(3) length(fa_list)]);
    parameter = parameter_list{k};
    output_parameter = output_parameter_list{k};
    for l = 1:length(fa_list)
        display(['k: ' num2str(k) ', l:' num2str(l)]);
        
        fa = fa_list{l};
        fa_folder = [parameter_folder 'fa-' fa '/'];
        
        parameter_map_nii = load_untouch_nii([fa_folder 'BrainSample-2_ses-03_' parameter '_fa-' fa '_20_directions_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
        parameter_map = (1./parameter_map_nii.img).* mask_erode;
        parameter_map(isnan(parameter_map)) = 0;
        
        parameter_map_all(:, :, :, l) = parameter_map ;
        parameter_map_nii.img = parameter_map;

        save_untouch_nii(parameter_map_nii, [fa_folder 'BrainSample-2_ses-03_' output_parameter '_fa-' fa '_20_directions_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    end
    
    parameter_map_mean = mean(parameter_map_all, 4);
    parameter_map_std = std(parameter_map_all, [], 4);
    
    parameter_map_nii.img = parameter_map_mean;
    save_untouch_nii(parameter_map_nii, [mean_folder 'BrainSample-2_ses-03_' output_parameter '_mean_20_directions_polyfit_cartesian_with_theta_mask_noise4_register.nii.gz']);
    parameter_map_nii.img = parameter_map_std;
    save_untouch_nii(parameter_map_nii, [std_folder 'BrainSample-2_ses-03_' output_parameter '_std_20_directions_polyfit_cartesian_with_theta_mask_noise4_register.nii.gz']);
end
