% Small images around samples

close all
clear

data_folder = '/project/3015069.01/derived/Porcine-2/ses-mri01/';
input_data_folder = [data_folder 'gre/highres/reference/']
output_data_folder = [data_folder 'concatenate_signals/'];
individual_folder = [output_data_folder 'individual_samples/'];

time = (linspace(2.44,50.73,12)*1e-3)';
edge = 10;

sample_mask_nii = load_untouch_nii('/project/3015069.01/derived/Porcine-2/ses-mri01/gre/highres/results/Porcine-2_gre_highres_ref_samples_final.nii.gz');
sample_mask = single(sample_mask_nii.img);
nb_sample = length(unique(sample_mask)) - 1;

v1 = load_nii_img_only([data_folder '/dwi/results/Porcine-2_dwi_orientation-6_NSA-9_lpca_eddy_dti_V1_2_ref.nii.gz']);

load([output_data_folder 'rotation_ref_2_orientations.mat']);

for k = 1:nb_sample
    small_sample_concatenate_signal_all_orientations{k} = [];
end
 
nb_orientation = 6;
length_one_orientation = 25;

mean_theta_registered = zeros(nb_sample, nb_orientation);

for k = 1:nb_orientation
    display(['nb orientation : ' num2str(k)])

    input_orientation_data_folder = [input_data_folder 'orientation-' num2str(k) '/'];
    
    input_magn = [input_orientation_data_folder 'Porcine-2_gre_highres_orientation-' num2str(k) '_part-mag_gradunwarp_2_ref.nii.gz'];
    input_phase = [input_orientation_data_folder 'sepia/Porcine-2_gre_highres_orientation-' num2str(k) '_unwrapped-phase.nii.gz'];
        
    magn = load_nii_img_only(input_magn);
    phase = load_nii_img_only(input_phase);
    
    nb_TE = size(magn, 4);
    
    for l = 1:nb_sample
        
        display(['nb sample : ' num2str(l)])
        
        sample =  single(sample_mask == l);
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
        sample_mask_nii.hdr.datatype = 16;
        sample_mask_nii.hdr.bitpix = 32;

        small_sample_magn = magn(min_x:max_x, min_y:max_y, min_z:max_z, :);
        small_sample_phase = phase(min_x:max_x, min_y:max_y, min_z:max_z, :);
        small_sample_v1 = v1(min_x:max_x, min_y:max_y, min_z:max_z, :);
 
        small_sample_real_poly = zeros([small_dims nb_TE]);
        small_sample_imag_poly = zeros([small_dims nb_TE]);
        
        small_v1_registered = zeros([small_dims 3]);
        
        nb_pixel = sum(small_sample_mask, 'all')
        orientations = zeros(nb_pixel, 3);
        
        it = 0; 
        for m = 1:small_dims(1)
            for n = 1:small_dims(2)
                for o = 1:small_dims(3)
                    if small_sample_mask(m, n, o)
                        
                        temp_magn = squeeze(small_sample_magn(m, n, o, :) / small_sample_magn(m, n, o, 1));
                        
                        small_v1_registered(m, n, o, :) = rotation(:, :, k) * squeeze(small_sample_v1(m, n, o, :));
                        
                        it = it + 1
                        orientations(it, :) = squeeze(small_v1_registered(m, n, o, :));
                        
                        temp_phase = squeeze(small_sample_phase(m, n, o, :));
                        poly_coeff = polyfit(time, temp_phase, 1);
                        phase_norm_poly = temp_phase - (time*poly_coeff(1) + poly_coeff(2));
                        
                        complex_norm_poly = temp_magn .* exp(1i*phase_norm_poly);
                                                
                        real_norm_poly = real(complex_norm_poly);
                        imag_norm_poly = imag(complex_norm_poly);
                        
                        small_sample_real_poly(m, n, o, :) = real_norm_poly;
                        small_sample_imag_poly(m, n, o, :) = imag_norm_poly;
                    end
                end
            end
        end
        theta_registered = acos(abs(small_v1_registered(:, :, :, 3))) .* small_sample_mask;

        avg_orientation = averageOrientationsLogEuclidean(orientations);       
        mean_theta_registered(l, k) = acos(abs(avg_orientation(3)));

        small_sample_concatenate_signal_all_orientations{l} = cat(4, small_sample_concatenate_signal_all_orientations{l}, theta_registered, small_sample_real_poly, small_sample_imag_poly);
    end   
end



h = figure;
for k = 1:7
    plot(mean_theta_registered(k,:), 'LineWidth', 2)
    hold on 
end
title('Theta angle per sample and orientation')
xlabel('orientation')
ylabel('theta')
set(gca, 'FontSize', 12)

origin_orient = [1; 0; 0];

for k = 1:6
proposed_orient(k,:) = rotation(:, :, k) * origin_orient;
end
proposed_theta = acos(abs(proposed_orient(:,3)))
plot(proposed_theta, 'LineWidth', 2)

leg = legend('1', '2', '3', '4', '5', '6', '7', 'proposed')
title(leg, 'sample')


for k = 1:nb_sample   
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
    
    sample_mask_nii.hdr.dime.dim(1) = 3;
    sample_mask_nii.hdr.dime.dim(2:4) = small_dims;
    sample_mask_nii.hdr.dime.datatype = 16;
    sample_mask_nii.hdr.dime.bitpix = 32;
    sample_mask_nii.img = single(small_sample_mask);  
    sample_mask_path = [individual_folder 'Porcine-2_ses-mri01_small_sample_' num2str(k) '_mask.nii.gz'];
    save_untouch_nii(sample_mask_nii, sample_mask_path)

    sample_mask_nii.hdr.dime.dim(1) = 4;
    sample_mask_nii.hdr.dime.dim(5) = nb_orientation * length_one_orientation;
    sample_mask_nii.img = single(small_sample_concatenate_signal_all_orientations{k});

    concatenate_sample_path = [individual_folder 'Porcine-2_ses-mri01_small_sample_' num2str(k) '_concatenate_signal_all_orientations_polyfit_cartesian_with_theta.nii.gz'];
    
    save_untouch_nii(sample_mask_nii, concatenate_sample_path)
end

