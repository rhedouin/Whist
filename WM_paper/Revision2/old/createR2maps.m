clear
close all

base_folder = ['/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/'];
concatenate_folder = [base_folder 'concatenate_signals_9_orientations/'];
parameter_folder = [concatenate_folder 'parameter_maps/BrainSample2LorentzinaCorrection/'];

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_mp2rage_orientation-4_t1_brain_2_BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1_brain_register_masked_mask_f05.nii.gz'));

fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60'};

for k = 1:length(fa_list)
    fa = fa_list{k};
    flip_angle_folder = [parameter_folder fa '/'];
    
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