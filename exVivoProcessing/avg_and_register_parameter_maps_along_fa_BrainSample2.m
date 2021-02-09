%test signal reconstruction
clear
close all

base_folder = ['/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/'];
ref_folder = [base_folder 'references/'];
concatenate_folder = [base_folder 'concatenate_signals_9_orientations/'];
parameter_folder = [concatenate_folder 'parameter_maps/BrainSample2LorentzinaCorrection/'];

parameter_list = {'FVF', 'gRatio', 'R2Myelin', 'R2IntraExtraAxonal', 'weight', 'xiMyelin'};

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_mask_coronal_73.nii.gz'));

fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60'};

for l = 1:length(fa_list)

    fa = fa_list{l};
    flip_angle_folder = [parameter_folder fa '/'];
    
    for k = 1:length(parameter_list)
        parameter = parameter_list{k};
        parameter_map = load_untouch_nii([flip_angle_folder 'BrainSample-2_ses-03_' parameter '_' fa '_polyfit_cartesian_with_theta_noise4.nii.gz']);
        parameter_map.img = permute(parameter_map.img, [1 3 2]);
        parameter_map.img = parameter_map.img(end:-1:1, end:-1:1, :);
        
        save_untouch_nii(parameter_map, [flip_angle_folder 'BrainSample-2_ses-03_' parameter '_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
        return;

    end
end

for k = 1:length(label_list)
    clear all_maps
    k
    it = 0;
    for fa = fa_list
        
        fa_folder = [base_folder 'fa-' num2str(fa) '/'];
        
        map_nii = load_untouch_nii([fa_folder prefix label_list{k} '_fa-' num2str(fa)  suffix '.nii.gz']);
        map = map_nii.img;
        
        it = it+1;
        all_maps(:,:,:,it) = map;   
    end
    
    mean_map = mean(all_maps, 4);
    std_map = std(all_maps,[], 4);
     
    map_nii.img = mean_map;
    mean_folder = [base_folder 'mean/'];
    mkdir(mean_folder);
    save_untouch_nii(map_nii, [mean_folder prefix label_list{k} '_mean' suffix '.nii.gz'])

    map_nii.img = std_map;
    std_folder = [base_folder 'std/'];
    mkdir(std_folder);
    save_untouch_nii(map_nii, [std_folder  prefix label_list{k} '_std' suffix  '.nii.gz'])
end

sepia




