% Create figure from brain for WM paper
clear
close all

base_folder = ['/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/'];
ref_folder = [base_folder 'references/'];
concatenate_folder = [base_folder 'concatenate_signals_3_orientations/'];

parameter_list = {'FVF', 'gRatio', 'R2Myelin', 'R2IntraExtraAxonal', 'weight', 'xiMyelin'};

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_mask_coronal_73.nii.gz'));

fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60'};

orientation_list = {'457', '157', '246', '245', '679', '269', '467', '579', '125', '234'};

fa = 'fa-35';

for l = 1:length(parameter_list)
    parameter = parameter_list{l}
    
    for k = 1:length(orientation_list)
        current_orientation = orientation_list{k}
        orientation_folder = [concatenate_folder current_orientation 'orientations/parameter_maps/BrainSample2LorentzinaCorrection/' fa '/'];
        
        parameter_map(:,:,:,k) = load_nii_img_only([orientation_folder 'BrainSample-2_ses-03_' parameter '_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
        
    end
    
    map_nii = load_untouch_nii([orientation_folder 'BrainSample-2_ses-03_R2IntraExtraAxonal_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    
    map_nii.img = mean(parameter_map,4);
    save_untouch_nii(map_nii, [concatenate_folder 'mean_' fa '/BrainSample-2_ses-03_' parameter '_mean_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    
    map_nii.img = std(parameter_map,[],4);
    save_untouch_nii(map_nii, [concatenate_folder 'std_' fa '/BrainSample-2_ses-03_' parameter '_std_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
end




