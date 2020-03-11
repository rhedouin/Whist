% This example simulates field perturbation and corresponding multi GRE
% signal from WM models
% Two single models are provided, one 2D model (by defaut) and one 3D model
% (that you can uncomment to run)
% You can load your own WM model create by createOneWMModelExample.m
 
clear
close all


% myelin (required: T2, xi, xa)
model_parameters.myelin.T2 = 15*1e-3;
% model_parameters.myelin.T1 = 500*1e-3;
model_parameters.myelin.proton_density= 0.5; 

model_parameters.myelin.xi = -0.1;  % myelin anisotropic susceptibility (ppm)
model_parameters.myelin.xa = -0.1;  % myelin isotropic susceptibility (ppm)

% intra axonal (required: T2)
model_parameters.intra_axonal.T2 = 50*1e-3;
% model_parameters.intra_axonal.T1 = 1.5;
model_parameters.intra_axonal.proton_density= 1; 
model_parameters.intra_axonal.xi= 0; 

% extra axonal (required: T2)
model_parameters.extra_axonal.T2 = 50*1e-3;
% model_parameters.extra_axonal.T1 = 1.5;
model_parameters.extra_axonal.proton_density= 1; 
model_parameters.extra_axonal.xi= 0; 

% TE (required)
% model_parameters.TE = (2:3:59)*1e-3;
model_parameters.TE = linspace(0.0001,0.08,100); 

% optional, needed to include T1 effect in signal weights
model_parameters.flip_angle = 20;
model_parameters.TR = 60*1e-3;
% 
model_parameters.include_proton_density = 1;
model_parameters.include_T1_effect = 0;

% main magnetic field strength in Tesla (required)
model_parameters.B0 = 3;
% magnetic field orientation (required)

base_model_folder = '/project/3015069.04/WM_Models/N400/';

model_list = [1, 3, 5, 7];
FVF_list = [10, 20, 30, 40, 50, 60, 70, 80];
theta_list = linspace(0, pi/2, 6);
nb_orientations = length(theta_list);

gRatio = 60;

for FVF = FVF_list
    for num = model_list
        model_path = [base_model_folder 'FVF' num2str(FVF) '_N400_train' num2str(num) '/FVF' num2str(FVF) '_gRatio' num2str(gRatio) '_N400_train' num2str(num)'];
        load(model_path)
        
        model_parameters.mask = mask;

        number_dims = ndims(mask);
        model_parameters.dims = size(mask);
        
        for k = 1:nb_orientations
            
            theta = theta_list(k);
            phi = 0;
            model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
            
            % Estimate the relative weights of each compartment
            model_parameters = computeCompartmentSignalWeight(model_parameters);
            
%             model_parameters.no_mask_tensor_map = 5;
            
            %%%%%%%%%% Simulate the field perturbation from the WM model and the multi GRE signals
            [signal_original, field] = simulateSignalFromModel(axon_collection, model_parameters);
            keyboard;
        end
    end
end










