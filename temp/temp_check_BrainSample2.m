% temp check BrainSample2
clear
close all

A = load_untouch_nii('local-field-lbv/BrainSample-2_ses-02_gre_orientation-1_fa-20_phase.nii.gz');
phase = A.img;

A = load_untouch_nii('local-field-lbv/BrainSample-2_ses-02_gre_orientation-1_fa-20_mask_qsm.nii.gz');
mask = A.img;

time = linspace(1.7, 35.25, 12)*1e-3;

norm_phase_poly = zeros(size(phase));
dims = size(phase);

for k = 1:dims(1)
    k
    for l = 1:dims(2)
        for m = 1:dims(3)
            if mask(k,l, m)
                unwrap_temp = squeeze(unwrap(phase(k,l,m,:)))';       
                unwrap_phase(k,l,m,:) = unwrap_temp;
                
                poly_coeff = polyfit(time, unwrap_temp, 1);
                temp_phase_poly = unwrap_temp - (time*poly_coeff(1) + poly_coeff(2));
                norm_phase_poly(k,l,m,:) = temp_phase_poly;
                
            end
        end
    end
end

A




