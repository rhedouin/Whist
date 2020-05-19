% Compute mean and std parameter map for WM paper

parameter_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals_3_orientations/';
fa = ['fa-20'];

mean_folder = [parameter_folder 'mean_fa-20/'];
std_folder = [parameter_folder 'std_fa-20/'];
mkdir(mean_folder)
mkdir(std_folder)

parameter_list = {'FVF', 'gRatio', 'T2Myelin', 'T2IntraExtraAxonal', 'weight', 'xiMyelin', 'R2Myelin', 'R2IntraExtraAxonal'};
% parameter_list = {'R2Myelin', 'R2IntraExtraAxonal'};

parameter_map_mean = 0;
parameter_map_example = load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals_9_orientations/parameter_maps/noise4/fa-05/BrainSample-2_ses-03_FVF_fa-05_20_directions_polyfit_cartesian_with_theta_noise4.nii.gz');

mask_map = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_gre_orientation-4_brain_mask_all_register.nii.gz'));
dims = size(mask_map);

mask_erode = mask_map;

list_orientations = [457, 157, 246, 245, 679, 269, 467, 579, 125, 234];


for k = 1:length(parameter_list)
    parameter = parameter_list{k};
    parameter_map_all = zeros([dims(1) dims(2) dims(3) length(list_orientations)]);
    
    for m = 1:length(list_orientations)
        current_orientation = list_orientations(m);
        orientation_folder = ['/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals_3_orientations/' num2str(current_orientation) 'orientations/parameter_maps/noise4/'];
        fa_folder = [orientation_folder fa '/'];
        
        parameter_map_nii = load_untouch_nii([fa_folder 'BrainSample-2_ses-03_' parameter '_' fa '_20_directions_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
        parameter_map = parameter_map_nii.img;
        parameter_map_all(:, :, :, m) = parameter_map ;      
    end
    
    parameter_map_mean = mean(parameter_map_all, 4);
    parameter_map_std = std(parameter_map_all, [], 4);

    parameter_map_nii.img = single(parameter_map_mean);
    save_untouch_nii(parameter_map_nii, [mean_folder 'BrainSample-2_ses-03_' parameter '_mean_' fa '_20_directions_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    parameter_map_nii.img = parameter_map_std;
    save_untouch_nii(parameter_map_nii, [std_folder 'BrainSample-2_ses-03_' parameter '_std_' fa '_20_directions_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    
end

