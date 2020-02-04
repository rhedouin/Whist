% Test signal
clear
close all

model_folder = '/project/3015069.04/WM_Models/N400/';
load([model_folder '/FVF10_N400_train1/AxonMap_FVF10_gRatio85_N400_train1.mat'])

mask = zeros(dims);
mask(round(1000/3):round(2*1000/3), round(1000/3):round(2*1000/3)) = 1;
imagesc(mask)

% parameters
B0 = 3;

xa = -0.1;
xi = -0.1;

T2out = 1000*1e-3;
T2myel = 15*1e-3;

weight = 0.5;
time = linspace(1.7,35.25,12)*1e-3;

field_direction = [0 1 0];

signal = simulateSignalFromModel(axon_collection, mask, xa, xi,  T2out, T2myel, weight, time, B0, field_direction);

magn = abs(signal);
phase = angle(signal);

signal_rep = repmat([0 magn phase(3:end)], [1 9])
plot([magn phase])
figure
plot(signal_rep)
hold on 
plot(data_signal_selection)
keyboard;

data_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals/';
data_path = [data_folder 'BrainSample-2_ses-03_all_orientations_fa-20_concatenate_signal_2_orientation-4.nii.gz'];

data_signal_values = load_nii_img_only(data_path);
temp_shape = size(data_signal_values);
data_shape = temp_shape(1:3);
data_signal_reshape = permute(data_signal_values,[4 1 2 3]);

CSF_list_pixel = {[53, 69, 65], [74, 70, 65], [54, 69, 61], [72, 67, 61]}

for  k = 1:4
data_signal_selection = data_signal_reshape(:, CSF_list_pixel{k}(1), CSF_list_pixel{k}(2), CSF_list_pixel{k}(3));

figure
plot(signal_rep)
hold on 
plot(data_signal_selection)
end


