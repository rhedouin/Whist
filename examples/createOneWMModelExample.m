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
your_folder = '/project/3015069.04/code/'; 
% location of the toolbox
axon_dictionary_path = [your_folder 'data/axonMediumDict.mat'];
% axon dictionary path
expected_FVF = 0.7;
% expected FVF of the final WM model

%%%%%%%%%% Options set with default values
%%%%%%%%%% White matter model 
options.number_of_axons = 400;
options.dims = [1000 1000]; 
options.mask = zeros(options.dims); 
options.mask(round(options.dims(1)/3):round(2*options.dims(1)/3), round(options.dims(2)/3):round(2*options.dims(2)/3)) = 1;

%%%%%%%%%% Axons packing 
options.max_FVF = 0.85;
options.max_iteration = 5000;
options.packing_speed = 0.5;

%%%%%%%%%%% Axons dispersion
options.dispersion_mode = 'spread'; 
options.tolerance = 0.01;

%%%%%%%%%%% Change g-ratio
options.expected_g_ratio = 0.6;

%%%%%%%%%%% Plot / save
options.plot_model = 1;
options.save_model = '/project/3015069.04/WM_Models/toto.mat';

createOne2DWMModel(axon_dictionary_path, expected_FVF, options);








