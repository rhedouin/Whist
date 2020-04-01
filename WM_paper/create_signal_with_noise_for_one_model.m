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
load([rotation_folder 'Porcine-2_lowres_rotation_ref_2_orientations.mat']);
load([rotation_folder '20_fiber_orientations_no_rotation.mat']);

%%%%%%%%%%% Set parameters
% myelin (required: T2, xi, xa)
model_parameters.myelin.T2 = 15*1e-3;
model_parameters.myelin.T1 = 300*1e-3;
model_parameters.myelin.proton_density= 0.5; 

model_parameters.myelin.xi = -0.1;  % myelin anisotropic susceptibility (ppm)
model_parameters.myelin.xa = -0.1;  % myelin isotropic susceptibility (ppm)

% intra axonal (required: T2)
model_parameters.intra_axonal.T2 = 50*1e-3;
model_parameters.intra_axonal.T1 = 1.5;
model_parameters.intra_axonal.proton_density= 1; 
model_parameters.intra_axonal.xi= 0; 

% extra axonal (required: T2)
model_parameters.extra_axonal.T2 = 50*1e-3;
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
model_parameters.TE = (1.8:3.2:56.2)*1e-3;
TE = model_parameters.TE;

model_parameters.TE = (1:85)*1e-3;
TE = model_parameters.TE ;

% optional, needed to include T1 effect in signal weights
model_parameters.flip_angle = 40;
model_parameters.TR = 60*1e-3;
 
model_parameters.include_proton_density = 1;
model_parameters.include_T1_effect = 1;

% Estimate the relative weights of each compartment
model_parameters = computeCompartmentSignalWeight(model_parameters);

FVF_round = 70;
gRatio_round = 65;

nb_models = 8;
nb_rotations = 6;
nb_sample_per_model = 125;
nb_TE = length(TE);

length_one_orientation = 2*nb_TE + 1;
length_total = length_one_orientation*nb_rotations;

noise_list = [0, 0.005, 0.01, 0.02, 0.04];

for noise = noise_list
    
    cube_nii = load_untouch_nii([dict_folder 'zeros_cube10.nii.gz']);
    M = zeros(nb_models*nb_sample_per_model, length_total);
    for k = 1:nb_models
        model_folder = [model_base_folder 'FVF' num2str(FVF_round) '_N400_train' num2str(k) '/'];
        load([model_folder 'FVF' num2str(FVF_round) '_gRatio' num2str(gRatio_round) '_N400_train' num2str(k) '.mat']);
     
        % mask (required)
        model_parameters.mask = mask;
        model_parameters.dims = size(mask);
        
        %%%%%%%%%% Simulate the field perturbation from the WM model and the multi GRE signals
        for l = 1:nb_rotations
            main_dir =  fiber_directions(1,:)';
            model_parameters.field_direction = rotations(:,:,l)*main_dir;
            theta = acos(abs(model_parameters.field_direction(3)));
            
            [signal_original, field] = simulateSignalFromModel(axon_collection, model_parameters);
            
            iteration = (k-1) * nb_sample_per_model;
            
            for m = 1:nb_sample_per_model
                
                iteration = iteration + 1
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
                
                M(iteration, (l-1)*length_one_orientation + 1 : l*length_one_orientation) = [theta_noisy  real_complex_polyfit imag_complex_polyfit];
            end
        end
    end
    
    M = reshape(M, [10 10 10 length_total]);
    
    cube_nii.img = single(M);
    if noise == 0.005
        save_untouch_nii(cube_nii, [dict_folder 'cube10_noise05.nii.gz'])
    else
        save_untouch_nii(cube_nii, [dict_folder 'cube10_noise' num2str(100*noise) '.nii.gz'])
    end
end







