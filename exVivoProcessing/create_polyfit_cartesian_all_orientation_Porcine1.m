% Ex vivo post processing 
% Concatenate signals
base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/';
dti_folder = [base_folder 'dwi/orientation-9/dti/'];
input_gre_folder = [base_folder 'gre_renaud/'];

fa_list = {'05','10','15','20','35','60'};
fa_list = {'05'};

lfa = length(fa_list);

for korient = 1:9
    orient = ['orientation-' num2str(korient)]
    input_orientation_folder = [input_gre_folder orient '/'];
    
    input_mask_path = [input_orientation_folder 'BrainSample-2_ses-03_gre_' orient '_brain_mask.nii.gz'];
    mask_nii = load_untouch_nii(input_mask_path);
    mask = single(mask_nii.img);
    dims = size(mask);
      
    input_theta_path = [input_orientation_folder 'BrainSample-2_ses-03_dwi_' orient '_theta.nii.gz'];
    theta_nii = load_untouch_nii(input_theta_path);
    theta = theta_nii.img;

    for kfa = 1:lfa
        tic()
        fa = ['fa-' fa_list{kfa}]
        input_fa_folder = [input_orientation_folder fa '/'];
            
        %% load normed magn and phase unwrapping
        input_normed_magn_path = [input_fa_folder 'BrainSample-2_ses-03_' orient '_' fa '_normed-magn-masked.nii.gz'];        
        magn_nii = load_untouch_nii(input_normed_magn_path);
        
        norm_magn = single(magn_nii.img);
                
        input_unwarped_phase_path = [input_fa_folder 'BrainSample-2_ses-03_' orient '_' fa '_unwrapped-phase-masked.nii.gz'];
        unwarped_phase_nii = load_untouch_nii(input_unwarped_phase_path);
        unwarped_phase = unwarped_phase_nii.img;
         
        total_dims = size(unwarped_phase);
        dims = total_dims(1:3);
        nb_TE = total_dims(4);
        
        %% magnitude and phase normalization
        norm_complex_poly = zeros([dims nb_TE]);
        norm_phase_poly = zeros([dims nb_TE]);

        norm_real_poly = zeros([dims nb_TE]);
        norm_imag_poly = zeros([dims nb_TE]);
        
        for k = 1:dims(1)
            k
            for l = 1:dims(2)
                for m = 1:dims(3)
                    if (mask(k,l,m) ~= 0)
                        temp_magn = squeeze(norm_magn(k,l,m,:))';                       
                        temp_phase = squeeze(unwarped_phase(k,l,m,:))';
                        
                        poly_coeff = polyfit(time, temp_phase, 1);
                        temp_phase_poly = temp_phase - (time*poly_coeff(1) + poly_coeff(2));
                        norm_phase_poly(k,l,m,:) = temp_phase_poly;
                         
                        [temp_real temp_imag] = pol2cart(temp_phase_poly, temp_magn);
                        
                        norm_real_poly(k,l,m,:) = temp_real;
                        norm_imag_poly(k,l,m,:) = temp_imag;
                        
                    end
                end
            end
        end
        
        concatenate_signal_cartesian_polyfit = cat(4, theta, norm_real_poly, norm_imag_poly);
        concatenate_signal_polar_polyfit = cat(4, theta, norm_magn, norm_phase_poly);
        
        %% Save all images
        output_normed_phase_poly_path = [input_fa_folder 'BrainSample-2_ses-03_' orient '_' fa '_normed-phase-polyfit-masked.nii.gz'];

        output_normed_real_poly_path = [input_fa_folder 'BrainSample-2_ses-03_' orient '_' fa '_normed-real-polyfit-masked.nii.gz'];
        output_normed_imag_poly_path = [input_fa_folder 'BrainSample-2_ses-03_' orient '_' fa '_normed-imag-polyfit-masked.nii.gz'];
        output_concatenate_signal_cartesian_polyfit_path = [input_fa_folder 'BrainSample-2_ses-03_' orient '_' fa '_concatenate_signal_cartesian_polyfit.nii.gz'];
        output_concatenate_signal_polar_polyfit_path = [input_fa_folder 'BrainSample-2_ses-03_' orient '_' fa '_concatenate_signal_polar_polyfit.nii.gz'];

        magn_nii.hdr.dime.datatype = 16;
        magn_nii.hdr.dime.bitpix = 32;      
        magn_nii.hdr.dime.dim(5) = 25;
        magn_nii.img = single(concatenate_signal_cartesian_polyfit);
        save_untouch_nii(magn_nii, output_concatenate_signal_cartesian_polyfit_path);

        magn_nii.hdr.dime.dim(5) = 25;
        magn_nii.img = single(concatenate_signal_polar_polyfit);
        save_untouch_nii(magn_nii, output_concatenate_signal_polar_polyfit_path);
        
        magn_nii.hdr.dime.dim(5) = 12;
        magn_nii.img = single(norm_phase_poly);
        save_untouch_nii(magn_nii, output_normed_phase_poly_path);
        
        magn_nii.hdr.dime.dim(5) = 12;
        magn_nii.img = single(norm_real_poly);
        save_untouch_nii(magn_nii, output_normed_real_poly_path);
        
        magn_nii.hdr.dime.dim(5) = 12;
        magn_nii.img = single(norm_imag_poly);
        save_untouch_nii(magn_nii, output_normed_imag_poly_path);
        toc()
    end   
end

    