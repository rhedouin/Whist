%test signal reconstruction
clear
close all

suffix = '_polyfit_cartesian_with_theta_noise4';

base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals_9_orientations/parameter_maps/BrainSample2LorentzinaCorrection/';

% loop all orientations
label_list = {'FVF', 'gRatio', 'xiMyelin', 'T2IntraExtraAxonal', 'T2Myelin', 'weight'};

fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60'};
for k = 1:length(label_list)
    clear all_maps
    k
    
    label =  label_list{k};
    it = 0;
    for l = 1:length(fa_list)
        fa = fa_list{l};
        fa_folder = [base_folder fa '/'];
        
        map_nii = load_untouch_nii([fa_folder 'BrainSample-2_ses-03_' label '_' fa '_polyfit_cartesian_with_theta_noise4.nii.gz']);
        map = map_nii.img;
        
        it = it+1;
        all_maps(:,:,:,it) = map;   
    end
    
    mean_map = mean(all_maps, 4);
    std_map = std(all_maps,[], 4);
     
    map_nii.img = mean_map;
    mean_folder = [base_folder 'mean/'];
    mkdir(mean_folder);    
    save_untouch_nii(map_nii, [mean_folder 'BrainSample-2_ses-03_' label '_mean_polyfit_cartesian_with_theta_noise4.nii.gz']);

    map_nii.img = std_map;
    std_folder = [base_folder 'std/'];
    mkdir(std_folder);
    save_untouch_nii(map_nii, [std_folder  'BrainSample-2_ses-03_' label '_std_polyfit_cartesian_with_theta_noise4.nii.gz']);

end






