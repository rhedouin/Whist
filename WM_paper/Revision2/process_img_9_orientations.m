clear;

base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/';
concatenate_folder = [base_folder 'concatenate_signals_9_orientations/'];
parameter_folder = [concatenate_folder 'parameter_maps/BrainSample2LorentzianCorrection/'];

cd(parameter_folder)

%%%%%%%%%%%%%%% Register images
input = [base_folder 'references/BrainSample-2_ses-03_mp2rage_orientation-4_t1_brain_2_BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1_brain'];
ref = [base_folder 'references/ref'];

ext = '_register.nii.gz';

fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60'};

parameter_list = {'FVF', 'gRatio', 'T2Myelin', 'T2IntraExtraAxonal', 'weight', 'xiMyelin'};

other = {};

for k = 1:length(fa_list)
     fa = fa_list{k};
     fa_folder = [parameter_folder fa '/'];
     
     for l = 1:length(parameter_list)
        parameter = parameter_list{l};
        other{end+1} = [fa_folder 'BrainSample-2_ses-03_' parameter '_' fa '_polyfit_cartesian_with_theta_noise4'];
     end
end
  
other{end+1} = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/references/BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1';

flirtAB2C(input, other, ref, ext)

%%%%%%%%%%%%%%% Create R2 maps 
mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_mp2rage_orientation-4_t1_brain_2_BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1_brain_register_masked_mask_f05.nii.gz'));

fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60'};

for k = 1:length(fa_list)
    fa = fa_list{k};
    fa_folder = [parameter_folder fa '/'];
    
    input_filename = [fa_folder 'BrainSample-2_ses-03_T2IntraExtraAxonal_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
    output_filename = [fa_folder 'BrainSample-2_ses-03_R2IntraExtraAxonal_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
    
    parameter_map = load_untouch_nii(input_filename);
    parameter_map.img = parameter_map.img.*(parameter_map.img > 0.00005).*mask;
    parameter_map.img = 1./parameter_map.img;
    parameter_map.img(isinf(parameter_map.img)) = 0;
    
    save_untouch_nii(parameter_map, output_filename)
    
    input_filename = [fa_folder 'BrainSample-2_ses-03_T2Myelin_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
    output_filename = [fa_folder 'BrainSample-2_ses-03_R2Myelin_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
    
    parameter_map = load_untouch_nii(input_filename);
    parameter_map.img = parameter_map.img.*(parameter_map.img > 0.00001).*mask;
    parameter_map.img = 1./parameter_map.img;
    parameter_map.img(isinf(parameter_map.img)) = 0;
    
    save_untouch_nii(parameter_map, output_filename)
    
end

%%%%%%%%%%%%%%% Create mean and std maps
parameter_list = {'FVF', 'gRatio', 'T2Myelin', 'T2IntraExtraAxonal', 'R2Myelin', 'R2IntraExtraAxonal', 'weight', 'xiMyelin'};

clear parameter_map

mkdir([parameter_folder 'mean'])
mkdir([parameter_folder 'std'])

for l = 1:length(parameter_list)
    parameter = parameter_list{l};

    for k = 1:length(fa_list)
        fa = fa_list{k};
        fa_folder = [parameter_folder fa '/'];
        
        parameter_map(:,:,:,k) = load_nii_img_only([fa_folder 'BrainSample-2_ses-03_' parameter '_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);        
    end
    
    map_nii = load_untouch_nii([fa_folder 'BrainSample-2_ses-03_R2IntraExtraAxonal_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    
    map_nii.img = mean(parameter_map,4);
    save_untouch_nii(map_nii, [parameter_folder 'mean/BrainSample-2_ses-03_' parameter '_mean_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    
    map_nii.img = std(parameter_map,[],4);
    save_untouch_nii(map_nii, [parameter_folder 'std/BrainSample-2_ses-03_' parameter '_std_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
end
