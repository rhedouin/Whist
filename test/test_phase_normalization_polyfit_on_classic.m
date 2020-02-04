% test_phase noise, background ...
clear
close all

model_folder = '/project/3015069.04/WM_Models/N400/';

figure
load([model_folder '/FVF70_N400_train6/AxonMap_FVF70_gRatio65_N400_train6.mat'])
%     load([model_folder '/FVF' num2str(k) '0_N400_train' num2str(k) '/AxonMap_FVF' num2str(k) '0_gRatio65_N400_train' num2str(k) '.mat'])

mask = zeros(dims);
mask(round(1000/3):round(2*1000/3), round(1000/3):round(2*1000/3)) = 1;

% parameters
B0 = 8;

xa = -0.1;
xi = -0.1;

T2out = 50*1e-3;
T2myel = 15*1e-3;

weight = 1;

% time = linspace(2.15,52.7,16)*1e-3;
time = linspace(2.15,35.7,12)*1e-3;

lTime = length(time);

figure

for theta = [0 pi/4 pi/2]
field_direction = [0 sin(theta) cos(theta)]

[signal_original, field] = simulateSignalFromModel(axon_collection, mask, xa, xi,  T2out, T2myel, weight, time, B0, field_direction);
zoomed_field = field(round(1000/3):round(2*1000/3), round(1000/3):round(2*1000/3));

phase_offset = 5;
freq_background = -10;


for k = 0:1:0
    
noise = k*0.005;
% noise = 0;
signal_noised = signal_original + noise*(randn(1,lTime) + 1i*randn(1,lTime));

signal_final = signal_noised .* exp(1i*(phase_offset + time * freq_background));
signal_final = signal_final / abs(signal_final(1));

real_original_signal = real(signal_final);
imag_original_signal = imag(signal_final);
magn_original_signal = abs(signal_final);
phase_original_signal = unwrap(angle(signal_final));

phase_norm_classic = phase_original_signal - phase_original_signal(1) - (phase_original_signal(2) - phase_original_signal(1))*(time - time(1))/(time(2) - time(1));
complex_norm_classic = magn_original_signal.*exp(1i*phase_norm_classic);
real_norm_classic = real(complex_norm_classic);
imag_norm_classic = imag(complex_norm_classic);

poly_coeff = polyfit(time, phase_original_signal, 1);
phase_norm_poly = phase_original_signal - (time*poly_coeff(1) + poly_coeff(2));
complex_norm_poly = magn_original_signal.*exp(1i*phase_norm_poly);
real_norm_poly = real(complex_norm_poly);
imag_norm_poly = imag(complex_norm_poly);
         

poly_coeff_bis = polyfit(time, phase_norm_classic, 1);
phase_norm_poly_bis = phase_norm_classic - (time*poly_coeff_bis(1) + poly_coeff_bis(2));
complex_norm_poly_bis = magn_original_signal.*exp(1i*phase_norm_poly_bis);
real_norm_poly_bis = real(complex_norm_poly_bis);
imag_norm_poly_bis = imag(complex_norm_poly_bis);

subplot(221)
hold on
plot([magn_original_signal phase_original_signal], 'LineWidth', 1.5)
title('original polar')

h = legend('0', 'pi/4', 'pi/2')

title(h, 'theta')
set(gca, 'FontSize', 15)

                
subplot(222)
hold on
plot([magn_original_signal phase_norm_classic], 'LineWidth', 1.5)
title('actual normed polar')
set(gca, 'FontSize', 15)

                     
subplot(223)
hold on
plot([magn_original_signal phase_norm_poly], 'LineWidth', 1.5)
title('polyfit normed polar')
set(gca, 'FontSize', 15)

         
subplot(224)
hold on
plot([magn_original_signal phase_norm_poly_bis], 'LineWidth', 1.5)
title('polyfit bis normed polar')
set(gca, 'FontSize', 15)

end
end

