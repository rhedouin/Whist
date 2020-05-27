% Transform Magn Phase in Polyfit cartesian
function [real_complex_polyfit, imag_complex_polyfit] =  TransformMagnPhase_2_PolyfitCartesian(magn, phase, mask)

TE = (1:size(phase, 4))';

dims = size(magn);
real_complex_polyfit = zeros(dims);
imag_complex_polyfit = zeros(dims);

for k = 1:dims(1)
    for l = 1:dims(2)
        for m = 1:dims(3)
            if mask(k, l, m)

                norm_magn = squeeze(magn(k, l, m, :)) / magn(k, l, m, 1);
                
                temp_phase = squeeze(phase(k, l, m, :));
                                
                poly_coeff = polyfit(TE, temp_phase, 1);
                norm_phase_polyfit = temp_phase - (TE*poly_coeff(1) + poly_coeff(2));
                                
                norm_complex_polyfit = norm_magn.*exp(1i*norm_phase_polyfit);
                real_complex_polyfit(k, l, m, :) = real(norm_complex_polyfit);
                imag_complex_polyfit(k, l, m, :) = imag(norm_complex_polyfit);
            end
        end
    end
end
end