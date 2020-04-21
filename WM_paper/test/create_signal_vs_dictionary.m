% Create nifti 10x10 from WM model
clear
close all
cd '/project/3015069.04/code/Whist';
your_folder = [pwd,'/']; 
% location of the toolbox
addpath(genpath(your_folder))

model_base_folder = '/project/3015069.04/WM_Models/N400/';
dict_folder = '/project/3015069.04/dictionaries/unique_value/';

rotation_folder = '/project/3015069.04/data/rotations/';
load([rotation_folder 'BrainSample2_rotations_ref_2_orientations.mat']);
load([rotation_folder '20_fiber_orientations_rotations.mat']);

%%%%%%%%%%% Set parameters
% myelin (required: T2, xi, xa)
model_parameters.myelin.T2 = 16*1e-3;
model_parameters.myelin.T1 = 300*1e-3;
model_parameters.myelin.proton_density= 0.5; 

model_parameters.myelin.xi = -0.1;  % myelin anisotropic susceptibility (ppm)
model_parameters.myelin.xa = -0.1;  % myelin isotropic susceptibility (ppm)

% intra axonal (required: T2)
model_parameters.intra_axonal.T2 = 60*1e-3;
model_parameters.intra_axonal.T1 = 1.5;
model_parameters.intra_axonal.proton_density= 1; 
model_parameters.intra_axonal.xi= 0; 

% extra axonal (required: T2)
model_parameters.extra_axonal.T2 = 60*1e-3;
model_parameters.extra_axonal.T1 = 1.5;
model_parameters.extra_axonal.proton_density= 1; 
model_parameters.extra_axonal.xi= 0; 

% main magnetic field strength in Tesla (required)
model_parameters.B0 = 3;
% magnetic field orientation (required)
theta = pi/2;
phi = 0;
model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];

% TE (required)
model_parameters.TE = (1.7:3.05:35.3)*1e-3;
TE = model_parameters.TE;

% optional, needed to include T1 effect in signal weights
model_parameters.flip_angle = 40;
model_parameters.TR = 60*1e-3;
 
model_parameters.include_proton_density = 0;
model_parameters.include_T1_effect = 0;

% Estimate the relative weights of each compartment
model_parameters = computeCompartmentSignalWeight(model_parameters);

model_parameters.myelin.weight = 1.5;
model_parameters.intra_axonal.weight = 1;
model_parameters.extra_axonal.weight = 1;

FVF_round = 70;
gRatio_round = 65;

nb_models = 8;
nb_rotations = 9;
nb_sample_per_model = 1;
nb_TE = length(TE);

length_one_orientation = 2*nb_TE + 1;
length_total = length_one_orientation*nb_rotations;

noise_list = [0, 0.005, 0.01, 0.02, 0.04];

dict_folder = '/project/3015069.04/dictionaries/multi_orientations/BrainSample2/';
dict_path = [dict_folder 'SignalWithNoise0_8rep_9rotations_12TE_BrainSample2_fix_xa_large_FVF_20_directions_polyfit_cartesian_with_theta.h5py'];

dict_signal_values = h5read(dict_path, '/SignalValues');
dict_signal_shape = size(dict_signal_values);
signal_length = dict_signal_shape(1);
dict_shape = dict_signal_shape(2:end);
dict_length = prod(dict_shape);

FVF_grid =    h5read(dict_path, '/FVFValues');
gRatio_grid = h5read(dict_path, '/gRatioValues');
xiMyelin_grid =     h5read(dict_path, '/xiMyelinValues');
T2IntraExtraAxonal_grid =    h5read(dict_path, '/T2IntraExtraAxonalValues');
T2Myelin_grid=    h5read(dict_path, '/T2MyelinValues');
weight_grid =    h5read(dict_path, '/weightValues');
direction_grid = h5read(dict_path, '/directionsValues');

dict_signal_reshape = reshape(dict_signal_values, signal_length, dict_length);
FVF_reshape = FVF_grid(:);
gRatio_reshape = gRatio_grid(:);
xiMyelin_reshape = xiMyelin_grid(:);
T2IntraExtraAxonal_reshape = T2IntraExtraAxonal_grid(:);
T2Myelin_reshape = T2Myelin_grid(:);
weight_reshape = weight_grid(:);

noise = 0;
k = 1;

model_folder = [model_base_folder 'FVF' num2str(FVF_round) '_N400_train' num2str(k) '/'];
load([model_folder 'FVF' num2str(FVF_round) '_gRatio' num2str(gRatio_round) '_N400_train' num2str(k) '.mat']);
     
% mask (required)
model_parameters.mask = mask;
model_parameters.dims = size(mask);
        
%%%%%%%%%% Simulate the field perturbation from the WM model and the multi GRE signals
clear data_signal
for l = 1:nb_rotations
    display(['nb rotations : ' num2str(l)])
    main_dir =  fiber_directions(1,:)';
    model_parameters.field_direction = rotations(:,:,l)*main_dir;
    theta = acos(abs(model_parameters.field_direction(3)));
    
    [signal_original, field] = simulateSignalFromModel(axon_collection, model_parameters);
    
    signal_noised = signal_original.total_normalized + noise*(randn(1,nb_TE) + 1i*randn(1,nb_TE));
    signal_noised = signal_noised / abs(signal_noised(1));
    
    norm_magn =  abs(signal_noised);
    phase = angle(signal_noised);
    
    poly_coeff = polyfit(TE, phase, 1);
    norm_phase_polyfit = phase - (TE*poly_coeff(1) + poly_coeff(2));
    
    norm_complex_polyfit= norm_magn.*exp(1i*norm_phase_polyfit);
    real_complex_polyfit = real(norm_complex_polyfit);
    imag_complex_polyfit = imag(norm_complex_polyfit);
    
    theta_noisy = theta + noise*randn(1);
    
    data_signal((l-1)*length_one_orientation + 1 : l*length_one_orientation,1) = [theta_noisy  real_complex_polyfit imag_complex_polyfit];
end

data_signal_replicate = repmat(data_signal,1,dict_length);
mse = sum((data_signal_replicate - dict_signal_reshape).^2, 1);

[mse_max, mse_sort] = sort(mse);

nb_map = 5;
for nb = 1:nb_map
    best_dict_match = dict_signal_reshape(:, mse_sort(nb));
    
    FVF(nb) = FVF_reshape(mse_sort(nb));
    gRatio(nb) = gRatio_reshape(mse_sort(nb));
    xiMyelin(nb) = xiMyelin_reshape(mse_sort(nb));
    T2IntraExtraAxonal(nb) = T2IntraExtraAxonal_reshape(mse_sort(nb));
    T2Myelin(nb) = T2Myelin_reshape(mse_sort(nb));
    weight(nb) = weight_reshape(mse_sort(nb));
    
    figure;
    plot(best_dict_match)
    hold on 
    plot(data_signal)

    legend('data', 'dict best match')

    title(['FVF : ' num2str(FVF(nb)) ', gRatio : ' num2str(gRatio(nb)) ', xi : ' num2str(xiMyelin(nb)) ...
         ', T2IntraExtraAxonal : ' num2str(T2IntraExtraAxonal(nb)) ', T2Myelin : ' num2str(T2Myelin(nb)) ', weight : ' num2str(weight(nb))])
    
end







