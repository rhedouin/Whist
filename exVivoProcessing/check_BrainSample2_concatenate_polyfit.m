% Test load and check polyfit Brain Sample 2
clear
close all
time = (1:12)';

lOrient = 23;
lOrient_polyfit = 25;

nb_orient = 9;

mask = load_nii_img_only('../masks/BrainSample-2_ses-03_gre_orientation-4_brain_mask.nii.gz');

fa_list = {'fa-60'};
for kFa = 1:length(fa_list)
    fa = fa_list{kFa};
    
    A = load_untouch_nii(['BrainSample-2_ses-03_all_orientations_' fa '_concatenate_signal_2_orientation-4.nii.gz']);
    signal = A.img;

    dims = size(signal);
    new_signal = zeros([dims(1) dims(2) dims(3) lOrient_polyfit * nb_orient]);
    
    for k = 1:dims(1)
        k
        for l = 1:dims(2)
            for m = 1:dims(3)
                if (mask(k,l,m)) == 1
                    new_temp_signal = [];
                    for orient = 1:nb_orient
                        temp_signal = squeeze(signal(k, l, m, :));

                        temp_theta = temp_signal((orient -1) * lOrient + 1);
                        temp_magn = temp_signal((orient -1) * lOrient + (2:13));
                        
                        temp_phase = temp_signal((orient -1) * lOrient + (14:23));
                        temp_phase = [0; 0; temp_phase];
                        
                        poly_coeff = polyfit(time, temp_phase, 1);
                        phase_norm_poly = temp_phase - (time*poly_coeff(1) + poly_coeff(2));
                        
                        complex_norm_poly = temp_magn .* exp(1i*phase_norm_poly);
                        
                        real_norm_poly = real(complex_norm_poly);
                        imag_norm_poly = imag(complex_norm_poly);
                        
                        new_temp_signal = [new_temp_signal; temp_theta; real_norm_poly; imag_norm_poly];
                        
                    end
                    new_signal(k, l, m, :) = new_temp_signal;
                end
            end
        end
    end
    A.img = single(new_signal);
    A.hdr.dime.dim(5) = lOrient_polyfit * nb_orient;
    
    save_untouch_nii(A, ['BrainSample-2_ses-03_all_orientations_' fa '_concatenate_signal_cartesian_polyfit_2_orientation-4.nii.gz']);
end
    
    
    
    
    
    
    
    
    
    
    
    
    