% Compute mean and std parameter map for WM paper

parameter_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals/parameter_maps/';
mean_folder = [parameter_folder 'mean/'];
std_folder = [parameter_folder 'std/'];

fa_list = {'05', '10', '15', '20', '35', '60'}

parameter_list = {'FVF', 'gRatio', 'T2Myelin', 'T2IntraExtraAxonal', 'weight', 'xiMyelin'};

parameter_map_mean = 0;
parameter_map_example = load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals/parameter_maps/fa-05_20_directions_polyfit_cartesian_with_theta/BrainSample-2_ses-03_FVF_fa-05_20_directions_polyfit_cartesian_with_theta.nii.gz');
dims = size(parameter_map_example);

mask_map = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_gre_orientation-4_brain_mask_int_without_external_csf.nii.gz'));
% se = strel('cube',3);
% mask_erode = imerode(mask_map,se);
mask_erode = mask_map;

for k = 1:length(parameter_list)
    parameter_map_all = zeros([dims(1) dims(2) dims(3) length(fa_list)]);
    parameter = parameter_list{k};
    for l = 1:length(fa_list)
        display(['k: ' num2str(k) ', l:' num2str(l)]);
        
        fa = fa_list{l};
        fa_folder = [parameter_folder 'fa-' fa '_20_directions_polyfit_cartesian_with_theta/'];
        
        parameter_map_nii = load_untouch_nii([fa_folder 'BrainSample-2_ses-03_' parameter '_fa-' fa '_20_directions_polyfit_cartesian_with_theta.nii.gz']);
        parameter_map = parameter_map_nii.img;
        
        parameter_map_all(:, :, :, l) = parameter_map .* mask_erode;
        
    end
    parameter_map_mean = mean(parameter_map_all, 4);
    parameter_map_std = std(parameter_map_all, [], 4);
    
    parameter_map_nii.img = parameter_map_mean;
    save_untouch_nii(parameter_map_nii, [mean_folder 'BrainSample-2_ses-03_' parameter '_mean_20_directions_polyfit_cartesian_with_theta_mask_without_external_csf.nii.gz']);
    parameter_map_nii.img = parameter_map_std;
    save_untouch_nii(parameter_map_nii, [std_folder 'BrainSample-2_ses-03_' parameter '_std_20_directions_polyfit_cartesian_with_theta_mask_without_external_csf.nii.gz']);
end
