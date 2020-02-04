% Small images around samples

close all
clear

data_folder = '/project/3015069.01/derived/Porcine-1/ses-mri02/concatenate_signals/';
input_data_folder = [data_folder 'individual_samples/'];

time = linspace(1.7,35.25,12)*1e-3;
edge = 10;

nb_sample = 7;

for k = 1:nb_sample
    small_sample_mask_path = ([input_data_folder 'Porcine-1_ses-mri02_small_sample_' num2str(k) '_mask.nii.gz']);
    small_sample_mask_nii = load_untouch_nii(small_sample_mask_path);
    small_sample_mask = small_sample_mask_nii.img;
    
    small_sample_path = ([input_data_folder 'Porcine-1_ses-mri02_small_sample_' num2str(k) '_concatenate_signal_all_orientations_without_theta_2_ref.nii.gz']);
    small_sample_nii = load_untouch_nii(small_sample_path);
    small_sample = single(small_sample_nii.img);

    small_dims = size(small_sample_mask);
    
    % Optional polyfit
    input_lVec = 23;
    output_lVec = 24;
    
    small_sample_polyfit = zeros([small_dims(1), small_dims(2), small_dims(3), 9*output_lVec]);
    for l = 1:small_dims(1)
        for m = 1:small_dims(2)
            for n = 1:small_dims(3)
                if small_sample_mask(l, m, n)
                    for orient = 1:9

                        temp_theta = small_sample(l, m, n, 1 + input_lVec*(orient - 1));
                        temp_magn = squeeze(small_sample(l, m, n, 2 + input_lVec*(orient - 1):13 + input_lVec*(orient - 1)))';

                        small_sample_polyfit(l, m, n, 1 + output_lVec*(orient - 1) :13 + output_lVec*(orient - 1)) = small_sample(l, m, n, 1 + input_lVec*(orient - 1):13 + input_lVec*(orient - 1));
                        
                        temp_phase = [0; 0; squeeze(small_sample(l, m, n, 14 + input_lVec*(orient - 1):23 + input_lVec*(orient - 1)))]';
                        
                        poly_coeff = polyfit(time, temp_phase, 1);
                        phase_norm_poly = temp_phase - (time*poly_coeff(1) + poly_coeff(2));
                        complex_norm_poly = temp_magn .* exp(1i*phase_norm_poly);
                        
                        real_norm_poly = real(complex_norm_poly);
                        imag_norm_poly = imag(complex_norm_poly);
                        
%                         small_sample_polyfit(l, m, n, 1 + output_lVec*(orient - 1)) = temp_theta;

                        small_sample_polyfit(l, m, n, 1 + output_lVec*(orient - 1) : 12 + output_lVec*(orient - 1)) = real_norm_poly;
                        small_sample_polyfit(l, m, n, 13 + output_lVec*(orient - 1) : 24 + output_lVec*(orient - 1)) = imag_norm_poly;
                        
                    end
                end
            end
        end
    end
    
    output_small_sample_path = ([input_data_folder 'Porcine-1_ses-mri02_small_sample_' num2str(k) '_concatenate_signal_polyfit_cartesian_all_orientations_without_theta_2_ref.nii.gz']);
    
    small_sample_nii.img = small_sample_polyfit;
    small_sample_nii.hdr.dime.dim(5) = 9*output_lVec;
    
    save_untouch_nii(small_sample_nii, output_small_sample_path);
    
    small_sample_mask_nii = 
    save_untouch_nii(small_sample_mask_nii, 
end



