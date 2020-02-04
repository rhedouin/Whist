% test_phase noise, background ...
clear
% close all

model_folder = '/project/3015069.04/WM_Models/N400/';

load([model_folder '/FVF70_N400_train6/AxonMap_FVF70_gRatio65_N400_train6.mat'])
%     load([model_folder '/FVF' num2str(k) '0_N400_train' num2str(k) '/AxonMap_FVF' num2str(k) '0_gRatio65_N400_train' num2str(k) '.mat'])

mask = zeros(dims);
mask(round(1000/3):round(2*1000/3), round(1000/3):round(2*1000/3)) = 1;

% parameters
B0 = 3;

xa = -0.1;
xi = -0.1;

T2out = 50*1e-3;
T2myel = 15*1e-3;

time = (2:3:59)*1e-3;

rho = 0.5;
TR = 65*1e-3;
T1out = 1.5; 
T1myel = 300 *1e-3;
flip_angle = 60;

weight = 1

figure
theta = pi/2;

field_direction = [0 sin(theta) cos(theta)];

phase_offset = 0;
freq_background = 0;

time_ms = time*1e3;
nbTE = length(time);

[signal_original, field] = simulateSignalFromModel(axon_collection, mask, xa, xi,  T2out, T2myel, weight, time, B0, field_direction);
signal_original = signal_original / abs(signal_original(1));

zoomed_field = field(round(1000/3):round(2*1000/3), round(1000/3):round(2*1000/3));

noise = 0.005;
signal_noised = signal_original + noise*(randn(1,nbTE) + 1i*randn(1,nbTE));

signal_final = signal_noised .* exp(1i*(phase_offset + time * freq_background));
signal_final = signal_final / abs(signal_final(1));

real_original_signal = real(signal_final);
imag_original_signal = imag(signal_final);
magn_original_signal = abs(signal_final);
phase_original_signal = unwrap(angle(signal_final));

poly_coeff = polyfit(time, phase_original_signal, 1);
phase_norm_poly = phase_original_signal - (time*poly_coeff(1) + poly_coeff(2));
complex_norm_poly = magn_original_signal.*exp(1i*phase_norm_poly);
real_norm_poly = real(complex_norm_poly);
imag_norm_poly = imag(complex_norm_poly);

figure;
hold on
plot(1:nbTE, real_norm_poly , 'LineWidth', 1.5, 'b')
plot(nbTE+3 : 2*nbTE+2, imag_norm_poly, 'LineWidth', 'b')
xline(nbTE + 1.5);
xticks([1, nbTE-1, nbTE+4, 2*nbTE+2])
xticklabels({time_ms(1),time_ms(end),time_ms(1),time_ms(end)})
xlabel('Echo time (ms)')

set(gca, 'FontSize', 15)







