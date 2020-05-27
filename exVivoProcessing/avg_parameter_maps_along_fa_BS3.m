%test signal reconstruction
clear
close all

suffix = '_polyfit_cartesian_with_theta_noise1';

base_folder = '/project/3015069.01/derived/BS-3/ses-mri03/gre/lowres/ref_space/concatenate_signals/parameter_maps/noise1/';

% loop all orientations
label_list = {'FVF', 'gRatio', 'xiMyelin', 'T2IntraExtraAxonal', 'T2Myelin', 'weight'};
prefix = 'BS-3_ses-mri03_acq-lowres_';
fa_list = [5, 10, 20, 30, 40, 50, 70];
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






