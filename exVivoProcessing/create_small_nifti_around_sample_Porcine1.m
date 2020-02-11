% Small images around samples

close all
clear

ses = 'ses-mri02';

data_folder = ['/project/3015069.01/derived/Porcine-1/' ses '/concatenate_signals/'];
output_data_folder = [data_folder 'individual_samples/'];
mkdir(output_data_folder)

sample_mask_nii = load_untouch_nii([data_folder 'Porcine-1_' ses '_anat_ref_samples.nii.gz']);
sample_mask = single(sample_mask_nii.img);
nb_sample = length(unique(sample_mask)) - 1;

time = linspace(1.7,35.25,12)*1e-3;
edge = 10;

theta_nii = load_nii_img_only(

for k = 1:nb_sample
    display(['nb sample : ' num2str(k)])

    sample =  single(sample_mask == k);
    dims = size(sample);
    
    list_sample = find(sample);
    [x_sample, y_sample, z_sample] = ind2sub(dims, list_sample);
    
    min_x = min(x_sample) - edge;
    min_y = min(y_sample) - edge;
    min_z = min(z_sample) - edge;
    
    max_x = max(x_sample) + edge;
    max_y = max(y_sample) + edge;
    max_z = max(z_sample) + edge;    
    
    small_sample_mask = sample(min_x:max_x, min_y:max_y, min_z:max_z);
    
    small_dims = size(small_sample_mask);
    
    sample_mask_nii.hdr.dime.dim(2:4) = size(small_sample_mask);
    sample_mask_nii.img = single(small_sample_mask);
    sample_mask_nii.hdr.dime.datatype = 16;
    sample_mask_nii.hdr.dime.bitpix = 32;
    small_sample_mask_path = ([output_data_folder 'Porcine-1_' ses '_small_sample_' num2str(k) '_mask.nii.gz']);
    save_untouch_nii(sample_mask_nii, small_sample_mask_path);
    
    input_lVec = 23;
    output_lVec = 25;
    nb_orient = 9;
    
    input_total_lVec = input_lVec*nb_orient;
    output_total_lVec = output_lVec*nb_orient;
    
    small_sample = [];
    small_sample_polyfit = zeros(small_dims(1), small_dims(2), small_dims(3), output_total_lVec);
    
    for orient = 1:nb_orient
        display(['nb orientation : ' num2str(orient)])
        
        concatenate_signal_nii = load_untouch_nii([data_folder 'Porcine-1_' ses '_orientation-' num2str(k) '_concatenate_signal_with_theta_2_ref.nii.gz']);
        small_signal_nii = concatenate_signal_nii;
        
        small_sample = cat(4, small_sample, concatenate_signal_nii.img(min_x:max_x, min_y:max_y, min_z:max_z, :));
        
        % Optional polyfit
        
        
        for l = 1:small_dims(1)
            for m = 1:small_dims(2)
                for n = 1:small_dims(3)
                    if small_sample_mask(l, m, n)
                        small_sample_polyfit(l, m, n,[1:13] + output_lVec*(orient - 1)) = squeeze(small_sample(l, m, n, [1:13] + input_lVec*(orient - 1)));
                        temp_phase = [0, 0, squeeze(small_sample(l, m, n, [14:23] + input_lVec*(orient - 1)))'];
                        
                        poly_coeff = polyfit(time, temp_phase, 1);
                        phase_norm_poly = temp_phase - (time*poly_coeff(1) + poly_coeff(2));
                        
                        small_sample_polyfit(l, m, n, [14:25] + output_lVec*(orient - 1)) = phase_norm_poly;
                        
                    end
                end
            end
        end
    end
    
    small_sample_mask_rep = repmat(small_sample_mask, [1 1 1 input_total_lVec]);
    
    small_signal_nii.hdr.dime.dim(2:5) = size(small_sample);
    small_signal_nii.img = small_sample .* small_sample_mask_rep;
    small_sample_path = ([output_data_folder 'Porcine-1_' ses '_small_sample_' num2str(k) '_concatenate_signal_all_orientations_with_theta_2_ref.nii.gz']);
    save_untouch_nii(small_signal_nii, small_sample_path);
    
    small_sample_mask_rep = repmat(small_sample_mask, [1 1 1 output_total_lVec]);
    
    small_signal_nii.img = small_sample_polyfit .* small_sample_mask_rep;;
    small_signal_nii.hdr.dime.dim(2:5) = size(small_sample_polyfit);
    small_sample_polyfit_path = ([output_data_folder 'Porcine-1_' ses '_small_sample_' num2str(k) '_concatenate_signal_polyfit_cartesian_all_orientations_with_theta_2_ref.nii.gz']);
    save_untouch_nii(small_signal_nii, small_sample_polyfit_path);
    
end













