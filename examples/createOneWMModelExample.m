%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example shows how to use the function createOne2DWMModel to create a
% 2D WM model with 3 compartments (intra-axonal, extra-axonal, myelin) using 
% realistic axon shapes.
% It follows 3 steps, packing axons, remove axons and change g-ratio (optional)
% See video complete_WM_Model_creation.mp4 for a demonstration.

% A dictionary of axon shapes is provided in this software. It was
% estimated from a spinal cord image of a dog using AxonSeg toolbox.
% You can furnish your own axon shapes dictionary if it respects the same format. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

%%%%%%%%%% Set parameters
%%%%%%%%%% Required

% your_folder = '/project/3015069.04/code/'; 
your_folder = [pwd,'/']; 
% location of the toolbox
axon_dictionary_path = [your_folder 'data/axonMediumDict.mat'];
addpath([your_folder 'createWMModel'])
% axon dictionary path

%%%%%%%%%% Options set with default values
%%%%%%%%%% White matter model 
model_params.number_of_axons = 100; % 
model_params.dims = [1000 1000]; 
% options.mask = zeros(options.dims); 
% options.mask(round(options.dims(1)/3):round(2*options.dims(1)/3), round(options.dims(2)/3):round(2*options.dims(2)/3)) = 1;

%%%%%%%%%% Axons packing 
model_params.max_FVF = 0.6;
model_params.max_iteration = 1000;
model_params.packing_speed = 2;

%%%%%%%%%%% Axons dispersion
model_params.expected_FVF = 0.1;
model_params.dispersion_mode = 'spread'; 
model_params.tolerance = 0.01;

%%%%%%%%%%% Change g-ratio
model_params.expected_g_ratio = 0.5;

%%%%%%%%%%% Plot / save
model_params.plot_model = 1;
model_params.save_model = [your_folder 'WMmodel/Toto.mat'];

[axon_collection, Model, ZoomedModel] = createOne2DWMModel(axon_dictionary_path, model_params);








