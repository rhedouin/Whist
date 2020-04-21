% Create figure from brain for WM paper
clear
close all

concatenate_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals/';
parameter_folder = [concatenate_folder 'parameter_maps/noise4/'];

parameter_list = {'FVF', 'gRatio', 'T2Myelin', 'T2IntraExtraAxonal', 'weight', 'xiMyelin'};

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_gre_orientation-4_brain_mask_all.nii.gz'));


fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60'}

for k = 1:length(fa_list)
    fa = fa_list{k};
    
    flip_angle_folder = [parameter_folder fa '/'];
    
    for k = 1:length(parameter_list)
        parameter = parameter_list{k};
        parameter_map.(parameter) = load_untouch_nii([flip_angle_folder 'BrainSample-2_ses-03_' parameter '_' fa '_20_directions_polyfit_cartesian_with_theta_noise4.nii.gz']);
        
        if ((k == 3) || (k == 4))
            parameter_map.(parameter) = parameter_map.(parameter)*1e3;
        end
        
        register_img_nii = parameter_map.(parameter);
        register_img_nii.img = permute(parameter_map.(parameter).img, [1 3 2]);
        register_img_nii.img = register_img_nii.img(:, end:-1:1, :);
   
        save_untouch_nii(register_img_nii, [flip_angle_folder 'BrainSample-2_ses-03_' parameter '_' fa '_20_directions_polyfit_cartesian_with_theta_noise4_reoriented.nii.gz'])
        keyboard;
    end
end



mean_folder = [parameter_folder 'mean/'];

for k = 1:length(parameter_list)
    parameter = parameter_list{k};
    parameter_map.(parameter) = load_untouch_nii([mean_folder 'BrainSample-2_ses-03_' parameter '_mean_20_directions_polyfit_cartesian_with_theta_noise4.nii.gz']);
    
    if ((k == 3) || (k == 4))
        parameter_map.(parameter) = parameter_map.(parameter)*1e3;
    end
    
    register_img_nii = parameter_map.(parameter);
    register_img_nii.img = permute(parameter_map.(parameter).img, [1 3 2]);
    register_img_nii.img = register_img_nii.img(:, end:-1:1, :);
    
    save_untouch_nii(register_img_nii, [mean_folder 'BrainSample-2_ses-03_' parameter '_mean_20_directions_polyfit_cartesian_with_theta_noise4_reoriented.nii.gz'])
end



std_folder = [parameter_folder 'std/'];

for k = 1:length(parameter_list)
    parameter = parameter_list{k};
    parameter_map.(parameter) = load_untouch_nii([std_folder 'BrainSample-2_ses-03_' parameter '_std_20_directions_polyfit_cartesian_with_theta_noise4.nii.gz']);
    
    if ((k == 3) || (k == 4))
        parameter_map.(parameter) = parameter_map.(parameter)*1e3;
    end
    
    register_img_nii = parameter_map.(parameter);
    register_img_nii.img = permute(parameter_map.(parameter).img, [1 3 2]);
    register_img_nii.img = register_img_nii.img(:, end:-1:1, :);
    
    save_untouch_nii(register_img_nii, [std_folder 'BrainSample-2_ses-03_' parameter '_std_20_directions_polyfit_cartesian_with_theta_noise4_reoriented.nii.gz'])
end
