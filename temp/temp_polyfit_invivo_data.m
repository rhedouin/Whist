% Separate concatenate signal

base_folder = '/project/3015069.04/data/InVivoProject/';
sub_list = 1:9


fa_list = {'05', '10', '20', '30', '40', '50', '70'};
lRot = 25;

nb_rotations = 9;

time = linspace(2.15,25.7, 12)'; 

for sub = sub_list
    sub_folder = [base_folder 'sub-0' num2str(sub) '/'];
    for kFa = 1:length(fa_list)
        fa = fa_list{kFa};
        fa_folder = [sub_folder 'fa-' fa '/'];

        display(['sub: ' num2str(sub) ', fa: ' fa]);
    
        signal_path = [fa_folder 'sub-0' num2str(sub) '_gre_fa-' fa '_concatenate_signal_theta_magn_phase.nii.gz'];
        signal_nii = load_untouch_nii(signal_path);
        signal = signal_nii.img;

        mask_path = [fa_folder 'sub-0' num2str(sub) '_gre_fa-' fa '_magn_unring_mask.nii.gz'];
        mask = load_nii_img_only(mask_path);
        
        dims = size(mask);
        
        signal_polyfit = zeros([dims(1) dims(2) dims(3) 25]);
        
        for k = 1:dims(1)
            k
            for l = 1:dims(2)
                for m = 1:dims(3)
                    if (mask(k,l,m) ~= 0)
                        temp_theta = signal(k,l,m,1);
                        temp_magn = squeeze(signal(k,l,m,2:13));                      
                        temp_phase = [0; 0; squeeze(signal(k,l,m,14:23))];
                        
                        poly_coeff = polyfit(time, temp_phase, 1);
                        temp_phase_poly = temp_phase - (time*poly_coeff(1) + poly_coeff(2));
                        norm_phase_poly = temp_phase_poly;
                         
                        [temp_real, temp_imag] = pol2cart(temp_phase_poly, temp_magn);
                        
                        signal_polyfit(k, l, m, :) = [temp_theta; temp_real; temp_imag];            
                        
                    end
                end
            end
        end
                
        signal_nii.hdr.dime.dim(5) = 25;
        signal_nii.img = single(signal_polyfit);
        
        output_signal_path = [fa_folder 'sub-0' num2str(sub) '_fa-' fa '_polyfit_cartesian_with_theta.nii.gz'];       
        save_untouch_nii(signal_nii, output_signal_path);
    end
end


