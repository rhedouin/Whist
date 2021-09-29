% This example simulates field perturbation and corresponding multi GRE
% signal from WM models
% Two single models are provided, one 2D model (by defaut) and one 3D model
% (that you can uncomment to run)
% You can load your own WM model create by createOneWMModelExample.m
 
clear
% close all

your_folder = [pwd,'/']; 
% location of the toolbox
addpath(genpath(your_folder))

%%%%%%%%%%%% Load a WM model with a single 2D axon
model_path = '/project/3015069.04/code/Whist/data/oneAxon2D.mat';

%%%%%%%%%%%% Load your WM model 
% model_path = '/project/3015069.04/code/Whist/WMmodel/MyWMmodel.mat';

%%%%%%%%%%%% Load a WM model with a single 3D axon
% model_path = '/project/3015069.04/code/Whist/data/oneAxon3D.mat';

load(model_path)

number_dims = ndims(mask);
model_parameters.dims = size(mask);

plot_model = 1;

% Create and plot your model
[model, zoomed_model, FVF, g_ratio] = createModelFromData(axon_collection, mask, plot_model);

%%%%%%%%%%% Set parameters
% mask (required)
model_parameters.mask = mask;

% myelin (required: T2, xi, xa)
model_parameters.myelin.weight= 0.5; 
model_parameters.myelin.T2 = 15*1e-3;
model_parameters.myelin.T1 = 500*1e-3;
model_parameters.myelin.proton_density= 0.5; 

model_parameters.myelin.xi = -0.1;  % myelin anisotropic susceptibility (ppm)
model_parameters.myelin.xa = -0.1;  % myelin isotropic susceptibility (ppm)

% intra axonal (required: T2)
model_parameters.intra_axonal.weight = 1;
model_parameters.intra_axonal.T2 = 50*1e-3;
model_parameters.intra_axonal.T1 = 1.5;
model_parameters.intra_axonal.proton_density= 1; 
model_parameters.intra_axonal.xi= 0; 

% extra axonal (required: T2)
model_parameters.extra_axonal.weight = 1;
model_parameters.extra_axonal.T2 = 50*1e-3;
model_parameters.extra_axonal.T1 = 1.5;
model_parameters.extra_axonal.proton_density= 1; 
model_parameters.extra_axonal.xi= 0; 

% main magnetic field strength in Tesla (required)
model_parameters.B0 = 3;
% magnetic field orientation (required)
theta = pi/2;
phi = 0;
model_parameters.theta = pi/2;
model_parameters.phi = 0;

model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];

% TE (required)
model_parameters.TE = (2:3:80)*1e-3;

% optional, needed to include T1 effect in signal weights
model_parameters.flip_angle = 20;
model_parameters.TR = 60*1e-3;
% 
model_parameters.include_proton_density = 1;
model_parameters.include_T1_effect = 1;
%%%%%%%%%%%

% Estimate the relative weights of each compartment
% model_parameters = computeCompartmentSignalWeight(model_parameters);

%%%%%%%%%% Simulate the field perturbation from the WM model and the multi GRE signals
[signal_original, field] = simulateSignalFromModel(axon_collection, model_parameters);

%%%%%%%%%% Plot frequency histogramm, field and signals
options.mask = mask;

createHistogramFieldPerturbation(model, field, options);

if number_dims == 2
    plot2DFieldAndSignal(field, signal_original.total, model_parameters.TE, model_parameters.field_direction)
elseif number_dims == 3
    plot3DFieldAndSignal(field, signal_original.total, model_parameters.TE, model_parameters.field_direction)
else
    error('dimension of the model should be 2 or 3');
end






