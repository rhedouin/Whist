% This example simulates field perturbation and corresponding multi GRE
% signal from WM models
% Two models are set provide, one 2D model (by defaut) and one 3D model
% (that you can simply uncomment to run)
 
clear
close all

%%%%%%%%%%%% Load a WM model with a single 2D axon
model_path = '/project/3015069.04/code/Whist/data/oneAxon2D.mat';

%%%%%%%%%%%% Load a WM model with a single 3D axon
% model_path = '/project/3015069.04/code/Whist/data/oneAxon3D.mat';

load(model_path)

number_dims = ndims(mask);
model_parameters.dims = size(mask);

% plot the WM model, the fiber volume fraction (FVF) is computed within the
% mask represented by the red rectangle
plot_model = 1;
[model, zoomed_model, FVF, g_ratio] = createModelFromData(oneAxon, mask, 1);

%%%%%%%%%%% Set parameters
% field strength in Tesla
model_parameters.B0 = 3;

% myelin 
model_parameters.myelin.T2 = 15*1e-3;
model_parameters.myelin.T1 = 500*1e-3;
model_parameters.myelin.proton_density= 0.5; 

model_parameters.myelin.xi = -0.1;  % myelin anisotropic susceptibility (ppm)
model_parameters.myelin.xa = -0.1;  % myelin isotropic susceptibility (ppm)

% intra axonal
model_parameters.intra_axonal.T2 = 50*1e-3;
model_parameters.intra_axonal.T1 = 1.5;
model_parameters.intra_axonal.proton_density= 1; 
model_parameters.intra_axonal.xi= 0; 

% extra axonal
model_parameters.extra_axonal.T2 = 50*1e-3;
model_parameters.extra_axonal.T1 = 1.5;
model_parameters.extra_axonal.proton_density= 1; 
model_parameters.extra_axonal.xi= 0; 

% TE 
model_parameters.TE = (2:3:59)*1e-3;

% Relative water signal received from myelin compare to intra/extra compartment (see
% computeCompartmentSignalWeight.m)
model_parameters.flip_angle = 20;
model_parameters.TR = 60*1e-3;

model_parameters.include_proton_density = 1;
model_parameters.include_T1_effect = 0;


model_parameters = computeCompartmentSignalWeight(model_parameters);

% magnetic field orientation
theta = pi/2;
phi = 0;
model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];

% mask
model_parameters.mask = mask;

%%%%%%%%%%%%% Simulate the field perturbation from the WM model and the multi GRE signals
[signal_original, field] = simulateSignalFromModel(oneAxon, model_parameters);

options.new_figure = 1;
options.mask = mask;

createHistogramFieldPerturbation(model, field, options);

magn_signal = abs(signal_original.total);
phase_signal = phase(signal_original.total);

%%%%%%%%%%%%% Plot field and signals
if number_dims == 2
    plot2DFieldAndSignal(field, signal_original.total, model_parameters.field_direction)
elseif number_dims == 3
    plot3DFieldAndSignal(field, signal_original.total, model_parameters.field_direction)
else
    error('dimension of the model should be 2 or 3');
end








