%test signal reconstruction
clear
close all

suffix = '_fix_xa';

base_folder = ['/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals/parameter_maps/SignalWithNoise05_8rep_with_theta_9orientations_BrainSample2_fix_xa_log_scale_with_normalization/'];

% loop all orientations
label_list = {'FVF', 'gRatio', 'xi', 'T2myel', 'T2out', 'weight'};
prefix = 'BrainSample-2_ses-03_';
fa_list = {'05', '10', '15', '20', '35', '60'};
for k = 1:length(label_list)
    clear all_maps
    k
    for kfa = 1:length(fa_list)
        fa = fa_list{kfa};
        
        fa_dir = [base_folder 'fa-' fa '/'];
        
        map_nii = load_untouch_nii([fa_dir prefix label_list{k} '_fa-' fa  suffix '.nii.gz']);
        map = map_nii.img;
        
        all_maps(:,:,:,kfa) = map;   
    end
    
    mean_map = mean(all_maps, 4);
    std_map = std(all_maps,[], 4);
     
    map_nii.img = mean_map;
    mean_folder = [base_folder 'mean/'];
    mkdir(mean_folder);
    save_untouch_nii(map_nii, [mean_folder prefix label_list{k} suffix '_mean.nii.gz'])

    map_nii.img = std_map;
    std_folder = [base_folder 'std/'];
    mkdir(std_folder);
    save_untouch_nii(map_nii, [std_folder  prefix label_list{k} suffix '_std.nii.gz'])
end






