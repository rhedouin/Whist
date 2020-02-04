%test signal reconstruction
clear
close all

ses = 'ses-mri02';
base_folder    = ['/project/3015069.01/derived/Porcine-1/' ses '/concatenate_signals/parameter_maps/single_orientation/'];

% loop all orientations
label_list = {'FVF', 'gRatio', 'T2out', 'xi', 'T2myel', 'weight', 'dispersion'};
prefix = ['Porcine-1_' ses];
for k = 1:length(label_list)
    
    for kori = 1:9
        orientation = ['orientation-' num2str(kori)]

        suffix = ['sample_mask_fix_xa_'];

        orient_dir = [base_folder  suffix orientation '/'];
        
        map_nii = load_untouch_nii([orient_dir prefix '_' label_list{k} '_' suffix orientation '.nii.gz']);
        map = map_nii.img;
        
        all_maps(:,:,:,kori) = map;   
    end
    
    mean_map = mean(all_maps, 4);
    std_map = std(all_maps,[], 4);
     
    map_nii.img = mean_map;
    save_untouch_nii(map_nii, [base_folder '/mean/'  prefix '_' label_list{k} '_' suffix 'mean.nii.gz'])

    map_nii.img = std_map;
    save_untouch_nii(map_nii, [base_folder '/std/'  prefix '_' label_list{k} '_' suffix 'std.nii.gz'])
end






