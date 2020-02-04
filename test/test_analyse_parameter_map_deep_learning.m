ses = 'ses-mri01';
base_folder = ['/project/3015069.01/derived/Porcine-1/' ses '/concatenate_signals/individual_samples/'];
data_folder = [base_folder 'parameter_maps/deep_learning_polyfit_cartesian_with_norm_with_reg/'];

parameter_list = {'FVF', 'gRatio', 'xi', 'T2out', 'T2myel', 'weight'};
kParameter = length(parameter_list);

nb_sample = 7;
for k = 1:nb_sample
    mask = load_nii_img_only([base_folder 'Porcine-1_ses-mri01_small_sample_' num2str(k)  '_mask.nii.gz']);

    fileID = fopen([data_folder 'Porcine-1_ses-mri01_small_sample_' num2str(k)  '_result_deep_learning.txt'], 'w');
    fprintf(fileID, '%s %s \n', '    ', 'Parameter avg'); 

    for l = 1:kParameter
        parameter = parameter_list{l};        
        parameter_map = load_nii_img_only([data_folder 'Porcine-1_ses-mri01_small_sample_' num2str(k) '_' parameter '_deep_learning_polyfit_cartesian_with_norm_with_reg.nii.gz']);
        
        parameter_mean = sum(parameter_map, 'all') / sum(mask, 'all')     
        
        fprintf(fileID, '%s %f \n', parameter, parameter_mean);
    end
    fclose(fileID);

end