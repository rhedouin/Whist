% Small images around samples

close all
clear

data_folder = '/project/3015069.01/derived/Porcine-1/ses-mri02/concatenate_signals/';
input_data_folder = [data_folder 'individual_samples/'];

nb_sample = 7;

for k = 1:nb_sample
    
    small_sample_path = ([input_data_folder 'Porcine-1_ses-mri02_small_sample_' num2str(k) '_concatenate_signal_polyfit_cartesian_all_orientations_with_theta_2_ref.nii.gz']);
    small_sample_nii = load_untouch_nii(small_sample_path);
    small_sample = single(small_sample_nii.img);

    % Optional polyfit
    remove_index = 1:25:225;
    correct_index = 1:225;
    correct_index(remove_index) = [];
    small_sample_without_theta = small_sample(:, :, :, correct_index);
    input_lVec = 25;
    output_lVec = 24;

    output_small_sample_path = ([input_data_folder 'Porcine-1_ses-mri02_small_sample_' num2str(k) '_concatenate_signal_polyfit_cartesian_all_orientations_without_theta_2_ref.nii.gz']);
    output_small_sample_nii = small_sample_nii;
    output_small_sample_nii.hdr.dime.dim(5) = 9*output_lVec;
    output_small_sample_nii.img = small_sample_without_theta;
    
    save_untouch_nii(output_small_sample_nii, output_small_sample_path);
end



