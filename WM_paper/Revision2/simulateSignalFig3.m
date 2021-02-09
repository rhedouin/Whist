% This example simulates field perturbation and corresponding multi GRE
% signal from WM models
% Two single models are provided, one 2D model (by defaut) and one 3D model
% (that you can uncomment to run)
% You can load your own WM model create by createOneWMModelExample.m
 
clear
close all

your_folder = [pwd,'/']; 
% location of the toolbox
addpath(genpath(your_folder))

%%%%%%%%%%%% Load a WM model with a single 2D axon
model_path = '/project/3015069.04/WM_Models/oneAxon/oneAxon.mat';

%%%%%%%%%%%% Load a WM model with a single 3D axon
% model_path = '/project/3015069.04/code/Whist/data/oneAxon3D.mat';

load(model_path)
axon_collection = oneAxon;

number_dims = ndims(mask);
model_parameters.dims = size(mask);

plot_model = 0;
% Create and plot your model
[model, zoomed_model, FVF, g_ratio] = createModelFromData(axon_collection, mask, plot_model);

%%%%%%%%%%% Set parameters
% mask (required)
model_parameters.mask = mask;

% myelin (required: T2, xi, xa)
model_parameters.myelin.T2 = 15*1e-3;
% model_parameters.myelin.T1 = 500*1e-3;
model_parameters.myelin.proton_density= 0.5; 
model_parameters.myelin.weight= 0.5; 

% intra axonal (required: T2)
model_parameters.intra_axonal.T2 = 50*1e-3;
% model_parameters.intra_axonal.T1 = 1.5;
model_parameters.intra_axonal.proton_density= 1; 
model_parameters.intra_axonal.xi= 0; 
model_parameters.intra_axonal.weight= 1; 

% extra axonal (required: T2)
model_parameters.extra_axonal.T2 = 50*1e-3;
% model_parameters.extra_axonal.T1 = 1.5;
model_parameters.extra_axonal.proton_density= 1; 
model_parameters.extra_axonal.xi= 0; 
model_parameters.extra_axonal.weight= 1; 

% main magnetic field strength in Tesla (required)
model_parameters.B0 = 3;
% magnetic field orientation (required)
theta = pi/2;
phi = 0;
model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];

% TE (required)
% model_parameters.TE = (2:3:59)*1e-3;
model_parameters.TE = linspace(0.0001,0.06,70); 

% Estimate the relative weights of each compartment
% model_parameters = computeCompartmentSignalWeight(model_parameters);

model_parameters.myelin.xi = -0.1;  % myelin anisotropic susceptibility (ppm)
model_parameters.myelin.xa = -0.1;  % myelin isotropic susceptibility (ppm)


%%%%%%%%%% Plot frequency histogramm, field and signals
for theta = [0 pi/4 pi/2]
    model_parameters.theta = theta;
    model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    model_parameters.current_dir = model_parameters.field_direction;
    
    %%%%%%%%%% Simulate the field perturbation from the WM model and the multi GRE signals
    [signal_original, field_original] = simulateSignalFromModel(axon_collection, model_parameters);
    [signal_corrected, field_corrected] = simulateSignalFromModelWithLorentzianCorrection(axon_collection, model_parameters);
    
    options.create_figure = 1;
    options.line_style = '-';
    options.mask = mask;
    options.edges = (-15:0.4:15);
    options.LineWidth = 2;
    options.xlim = [-15 15];
    createHistogramFieldPerturbation(model, field_corrected, options);
    
    options.create_figure = 0;
    options.line_style = '--';
    set(gca, 'FontSize', 20, 'FontWeight','bold' )
    createHistogramFieldPerturbation(model, field_original, options);

end





