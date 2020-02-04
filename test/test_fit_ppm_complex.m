% test fit_ppm complex
clear 
close all

model_folder = '/project/3015069.04/WM_Models/N400/';

load([model_folder '/FVF70_N400_train6/AxonMap_FVF70_gRatio65_N400_train6.mat'])
%     load([model_folder '/FVF' num2str(k) '0_N400_train' num2str(k) '/AxonMap_FVF' num2str(k) '0_gRatio65_N400_train' num2str(k) '.mat'])

mask = zeros(dims);
mask(round(1000/3):round(2*1000/3), round(1000/3):round(2*1000/3)) = 1;

% parameters
B0 = 7;

xa = -0.1;
xi = -0.1;

T2out = 50*1e-3;
T2myel = 15*1e-3;

weight = 1;

% time = linspace(2.15,52.7,16)*1e-3;
% time = linspace(2.15,35.7,12)*1e-3;
time = linspace(1,12,12);

lTime = length(time);

field_direction = [0, 1, 0];
[signal_original, field] = simulateSignalFromModel(axon_collection, mask, xa, xi,  T2out, T2myel, weight, time, B0, field_direction);
zoomed_field = field(round(1000/3):round(2*1000/3), round(1000/3):round(2*1000/3));

imagesc(zoomed_field)

phase_offset = 0;
freq_background = 0;
    
figure
for k = 1:6
noise = (k-1)*0.0025;
signal_noised = signal_original + noise*(randn(1,lTime) + 1i*randn(1,lTime));

signal_final = signal_noised .* exp(1i*(phase_offset + time * freq_background));
signal_final = signal_final / abs(signal_final(1));

magn_original_signal = abs(signal_final);
phase_original_signal = unwrap(angle(signal_final));
complex_original_signal = magn_original_signal .* exp(1i*phase_original_signal);
rep_complex_original_signal = reshape(complex_original_signal, [1 1 1 lTime]);

poly_coeff = polyfit(time, phase_original_signal, 1)
phase_norm_poly = phase_original_signal - (time*poly_coeff(1) + poly_coeff(2));
complex_norm_poly = magn_original_signal.*exp(1i*phase_norm_poly);
% rep_complex_norm_poly = reshape(complex_norm_poly, [10 10 10 lTime]);
rep_complex_norm_poly = permute(repmat(squeeze(complex_norm_poly)', [1 5 5 5]), [2 3 4 1]);

% [~, ~, ~, phase_fit_ppm_1] = Fit_ppm_complex(rep_complex_norm_poly)
[p1, dp1, relres, p0] = Fit_ppm_complex(rep_complex_original_signal)
keyboard;
phase_fit_ppm = phase_original_signal - (time*poly_coeff(1) + poly_coeff(2));


keyboard;
subplot(2,3,k)
plot(phase_original_signal)
hold on

plot(phase_norm_poly)
plot(-phase_fit_ppm_1)
plot(-phase_fit_ppm_2)

legend('original phase', 'phase norm poly', 'phase fit ppm on poly', 'phase fit ppm on original')
keyboard;
end





