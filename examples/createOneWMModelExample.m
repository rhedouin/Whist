% createOneWMModel example

% set parameters
% required
axon_dictionary_path = '/project/3015069.04/histology/axon_dictionary/axonLargeDict.mat';
expected_FVF = 0.7;

% optional
% White matter model 
options.number_of_axons = 400;
options.side = 1000;

% Axons packing 
options.max_FVF = 0.85;
options.max_iteration = 5000;
options.packing_speed = 0.5;

% Axons dispersion
options.dispersion_mode = 'spread'; 
options.tolerance = 0.1;

% Change g-ratio
options.expected_g_ratio = 0.6;

% Plot / save
options.plot = 1;
options.save_model = '/project/3015069.04/WM_Models/toto.mat';


createOne2DWMModel(axon_dictionary_path, expected_FVF, options)








