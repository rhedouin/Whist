% Compute mean and std parameter map for WM paper

parameter_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals_9_orientations/parameter_maps/noise4/';
mean_folder = [parameter_folder 'mean/'];
std_folder = [parameter_folder 'std/'];

fa_list = {'05', '10', '15', '20', '35', '60'}

parameter_list = {'FVF', 'gRatio', 'T2Myelin', 'T2IntraExtraAxonal', 'weight', 'xiMyelin', 'T2Myelin', 'T2IntraExtraAxonal', 'R2Myelin', 'R2IntraExtraAxonal'};
% parameter_list = {'R2Myelin', 'R2IntraExtraAxonal'};

parameter_map_mean = 0;
parameter_map_example = load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals_9_orientations/parameter_maps/noise4/fa-05/BrainSample-2_ses-03_FVF_fa-05_20_directions_polyfit_cartesian_with_theta_noise4.nii.gz');

mask_map = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_gre_orientation-4_brain_mask.nii.gz'));
dims = size(mask_map);

mask_erode = mask_map;

list_orientations = [457, 157, 246, 245, 679, 269, 467, 579, 125, 234];
fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60'}
fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35'}

for m = 1:length(list_orientations)
    current_orientation = list_orientations(m)
    orientation_folder = ['/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals_3_orientations/' num2str(current_orientation) 'orientations/parameter_maps/noise4/'];
    mean_folder = [orientation_folder 'mean/'];
    std_folder = [orientation_folder 'std/'];
    
    mkdir(mean_folder)
    mkdir(std_folder)
        
    for k = 1:length(parameter_list)
        parameter_map_all = zeros([dims(1) dims(2) dims(3) length(fa_list)]);
        parameter = parameter_list{k};

        for l = 1:length(fa_list)
            display(['k: ' num2str(k) ', l:' num2str(l)]);
            
            fa = fa_list{l};
            fa_folder = [parameter_folder fa '/'];
            
            parameter_map_nii = load_untouch_nii([fa_folder 'BrainSample-2_ses-03_' parameter '_' fa '_20_directions_polyfit_cartesian_with_theta_noise4.nii.gz']);
            parameter_map = parameter_map_nii.img;

            parameter_map(isnan(parameter_map)) = 0;
            
            parameter_map_all(:, :, :, l) = parameter_map ;
            parameter_map_nii.img = parameter_map;
            
        end
        
        parameter_map_mean = mean(parameter_map_all, 4);
        parameter_map_std = std(parameter_map_all, [], 4);
    
        parameter_map_nii.img = parameter_map_mean;
        save_untouch_nii(parameter_map_nii, [mean_folder 'BrainSample-2_ses-03_' parameter '_mean_20_directions_polyfit_cartesian_with_theta_noise4.nii.gz']);
        parameter_map_nii.img = parameter_map_std;
        save_untouch_nii(parameter_map_nii, [std_folder 'BrainSample-2_ses-03_' parameter '_std_20_directions_polyfit_cartesian_with_theta_noise4.nii.gz']);
   
    end
end
