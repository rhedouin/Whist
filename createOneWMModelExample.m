%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example shows how to use the function createOne2DWMModel to create a
% 2D WM model with 3 compartments (intra-axonal, extra-axonal, myelin) using 
% realistic axon shapes.
% It follows 3 steps, packing axons, remove axons and change g-ratio (optional)
% See video complete_WM_Model_creation.mp4 for a demonstration.

% A dictionary of axon shapes is provided in this software. It was
% estimated from a spinal cord image of a dog using AxonSeg toolbox.
% You can furnish your own axon shapes dictionary if it respects the same format 
% (see createOne2DWMModel.m, axon_collection) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

%%%%%%%%%% Set parameters
%%%%%%%%%% Required

your_folder = '/project/3015069.04/code/Whist/'; 
% location of the toolbox
axon_dictionary_path = [your_folder 'data/axonMediumDict.mat'];
addpath([your_folder 'createWMModel'])
% axon dictionary path

%%%%%%%%%% Options 
%%%%%%%%%% White matter model 
model_params.number_of_axons = 400; % 
% 
%%%%%%%%%% Axons packing 
model_params.max_FVF = 0.8;
model_params.max_iteration = 5000;
model_params.packing_speed = 1;
% 
%%%%%%%%%%% Axons dispersion
model_params.expected_FVF = 0.5405;
model_params.dispersion_mode = 'spread'; 
model_params.tolerance = 0.01;

%%%%%%%%%%% Change g-ratio
model_params.expected_g_ratio = 0.6723;

%%%%%%%%%%% Plot / save
model_params.plot_model = 0;

for k = 3:10
    model_params.save_model = [your_folder 'WMmodel/2DModel_FVF054_gRatio_067_v' num2str(k) '.mat'];

    [axon_collection, Model, ZoomedModel] = createOne2DWMModel(axon_dictionary_path, model_params);
end






