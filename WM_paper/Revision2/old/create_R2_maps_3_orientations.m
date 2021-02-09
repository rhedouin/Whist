% Create figure from brain for WM paper
clear
close all

base_folder = ['/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/'];
ref_folder = [base_folder 'references/'];
parameter_folder = [base_folder 'concatenate_signals_3_orientations/'];

parameter_list = {'FVF', 'gRatio', 'T2Myelin', 'T2IntraExtraAxonal', 'weight', 'xiMyelin'};

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_mask_coronal_73.nii.gz'));

fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60'};
orientation_list = {'457', '157', '246', '245', '679', '269', '467', '579', '125', '234'};
 
for k = 1:length(orientation_list)
        
    current_orientation = orientation_list{k}
    orientation_folder = [parameter_folder current_orientation 'orientations/parameter_maps/BrainSample2LorentzinaCorrection/'];

    for l = 1:length(fa_list)
        
        fa = fa_list{l}
        flip_angle_folder = [orientation_folder fa '/'];
        
        input_filename = [flip_angle_folder 'BrainSample-2_ses-03_T2IntraExtraAxonal_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
        output_filename = [flip_angle_folder 'BrainSample-2_ses-03_R2IntraExtraAxonal_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
        
        parameter_map = load_untouch_nii(input_filename);
        parameter_map.img = parameter_map.img.*(parameter_map.img > 0.00005).*mask;
        parameter_map.img = 1./parameter_map.img;
        parameter_map.img(isinf(parameter_map.img)) = 0;
        
        save_untouch_nii(parameter_map, output_filename)
        
        input_filename = [flip_angle_folder 'BrainSample-2_ses-03_T2Myelin_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
        output_filename = [flip_angle_folder 'BrainSample-2_ses-03_R2Myelin_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
        
        parameter_map = load_untouch_nii(input_filename);
        parameter_map.img = parameter_map.img.*(parameter_map.img > 0.00001).*mask;
        parameter_map.img = 1./parameter_map.img;
        parameter_map.img(isinf(parameter_map.img)) = 0;

        save_untouch_nii(parameter_map, output_filename)        
    end
end

