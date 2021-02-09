clear;


base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/';
concatenate_folder = [base_folder 'concatenate_signals_3_orientations/'];
% experiment_name = 'BrainSample2LorentzianCorrection';
experiment_name = 'BrainSample2LorentzianCorrection';

fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60'};
parameter_list = {'FVF', 'gRatio', 'R2Myelin', 'R2IntraExtraAxonal', 'weight', 'xiMyelin'};
list_orientations = {'457', '157', '246', '245', '679', '269', '467', '579', '125', '234'};

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_mp2rage_orientation-4_t1_brain_2_BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1_brain_register_masked_mask_f05.nii.gz'));

%%%%%%%%%%%%%%% Register images
% input = [base_folder 'references/BrainSample-2_ses-03_mp2rage_orientation-4_t1_brain_2_BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1_brain'];
% ref = [base_folder 'references/ref'];
% 
% ext = '_register.nii.gz';
% 
% other = {};
% 
% for m = 1 : length(list_orientations)
%     current_orientation = list_orientations{m};
%     orientation_folder = [concatenate_folder current_orientation 'orientations/parameter_maps/' experiment_name '/'];
%     
%     for k = 1:length(fa_list)
%         fa = fa_list{k};
%         fa_folder = [orientation_folder fa '/'];
%         
%         for l = 1:length(parameter_list)
%             parameter = parameter_list{l};
%             other{end+1} = [fa_folder 'BrainSample-2_ses-03_' parameter '_' fa '_polyfit_cartesian_with_theta_noise4'];
%         end
%     end
% end
% 
% flirtAB2C(input, other, ref, ext)

% %%%%%%%%%%%%%% Create R2 maps 
% 
% for m = 1 : length(list_orientations)
%     current_orientation = list_orientations{m};
%     orientation_folder = [concatenate_folder current_orientation 'orientations/parameter_maps/' experiment_name '/'];
%     
%     for k = 1:length(fa_list)
%         fa = fa_list{k};
%         fa_folder = [orientation_folder fa '/'];
%         
%         input_filename = [fa_folder 'BrainSample-2_ses-03_T2IntraExtraAxonal_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
%         output_filename = [fa_folder 'BrainSample-2_ses-03_R2IntraExtraAxonal_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
%         
%         parameter_map = load_untouch_nii(input_filename);
%         parameter_map.img = parameter_map.img.*(parameter_map.img > 0.00005).*mask;
%         parameter_map.img = 1./parameter_map.img;
%         parameter_map.img(isinf(parameter_map.img)) = 0;
%         
%         save_untouch_nii(parameter_map, output_filename)
%         
%         input_filename = [fa_folder 'BrainSample-2_ses-03_T2Myelin_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
%         output_filename = [fa_folder 'BrainSample-2_ses-03_R2Myelin_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz'];
%         
%         parameter_map = load_untouch_nii(input_filename);
%         parameter_map.img = parameter_map.img.*(parameter_map.img > 0.00001).*mask;
%         parameter_map.img = 1./parameter_map.img;
%         parameter_map.img(isinf(parameter_map.img)) = 0;
%         save_untouch_nii(parameter_map, output_filename)
%     end
% end

%%%%%%%%%%%%%%% Create mean and std maps
parameter_list = {'FVF', 'gRatio', 'R2Myelin', 'R2IntraExtraAxonal', 'weight', 'xiMyelin'};

clear parameter_map

for k = 1:length(fa_list)
    fa = fa_list{k};
    mean_fa_folder = [concatenate_folder 'mean_' fa '/' experiment_name];
    std_fa_folder = [concatenate_folder 'std_' fa '/' experiment_name];
    
    mkdir(mean_fa_folder)
    mkdir(std_fa_folder)
    
    for l = 1:length(parameter_list)
        parameter = parameter_list{l};
        
        clear parameter_map
        
        for m = 1 : length(list_orientations)
            current_orientation = list_orientations{m};
            orientation_folder = [concatenate_folder current_orientation 'orientations/parameter_maps/' experiment_name '/'];
            fa_folder = [orientation_folder fa '/'];
            
            parameter_map(:,:,:,m) = load_nii_img_only([fa_folder 'BrainSample-2_ses-03_' parameter '_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
        end
        
        map_nii = load_untouch_nii([fa_folder 'BrainSample-2_ses-03_R2IntraExtraAxonal_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
        
        map_nii.img = mean(parameter_map,4);
        save_untouch_nii(map_nii, [mean_fa_folder '/BrainSample-2_ses-03_' parameter '_mean_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
        
        map_nii.img = std(parameter_map,[],4);
        save_untouch_nii(map_nii, [std_fa_folder '/BrainSample-2_ses-03_' parameter '_std_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    end
end
