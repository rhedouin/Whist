function [axon_collection, Model, ZoomedModel] = createOne2DWMModel(axon_dictionary_path, expected_FVF, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Inputs required
%  axon_dictionary_path   : input axon dictionary path (provided in the
%  data folder)
%  expected_FVF : FVF of the final WM model
 
% %%%%%%%%%% Options
% %%%%%%%%%% White matter model 
% options.number_of_axons = 400;
% options.dims = [1000 1000]; 
% % dims of the original axon grid (estimated from the number of axons as
% % well as the axon size if not provided)
% options.mask = zeros(options.dims); 
% options.mask(round(options.dims(1)/3):round(2*options.dims(1)/3), round(options.dims(2)/3):round(2*options.dims(2)/3)) = 1;
% % the mask is the area where the fiber volume fraction (FVF) will be
% % computed. It represents the actual final WM model after axon packing
% 
% %%%%%%%%%% Axons packing 
% options.max_FVF = 0.85;
% % FVF value where the axon packing is stopped (it cannot be much higher than
% % 0.85 )
% options.max_iteration = 5000;
% % Max number of iteration of the axon packing (stop the axon packing if it
% % cannot reach the max_FVF)
% options.packing_speed = 0.5;
% % The packing speed weigths the attraction/repulsion of the axons. A higher
% % value accelerates the packing process but can create axons overlap
% 
% %%%%%%%%%%% Axons dispersion
% options.dispersion_mode = 'spread'; 
% % Dispersion mode can be remove or spread
% options.tolerance = 0.001;
% % Tolerance between expected FVF and actual FVF of the model
% 
% %%%%%%%%%%% Change g-ratio
% options.expected_g_ratio = 0.6;
% % Change the myelin thickness to reach an expected g-ratio
% 
% %%%%%%%%%%% Plot / save model
% options.plot_model = 1;
% options.save_model = '/project/3015069.04/WM_Models/toto.mat';
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default parameters  
if ~exist(options)
    options.null = 0;
end

if ~isfield(options, 'max_FVF')
    options.max_FVF = 0.85;
end

if ~isfield(options, 'max_iteration')
    options.max_iteration = 5000;
end

if ~isfield(options, 'packing_speed')
    options.packing_speed = 0.5;
end

if ~isfield(options, 'dispersion_mode')
    options.dispersion_mode = 'spread';
end
 
if ~isfield(options, 'number_of_axons')
    options.number_of_axons = 400;
end

if ~isfield(options, 'plot_model')
    options.plot_model = 1;
end

if ~isfield(options, 'tolerance')
    options.tolerance = 0.1;
end

% Load dictionary
disp('load axon dictionary ...')
load(axon_dictionary_path);
disp('done')

disp(['randomly select ' num2str(options.number_of_axons) ' axon shapes ...']);
list_axons = randi(length(axonDico), options.number_of_axons, 1);
original_axon_collection = axonDico(list_axons);
for k = 1:length(original_axon_collection)
    original_axon_collection(k).data = double(original_axon_collection(k).data);
end
disp('done')

% Setup axons grid
disp('setup axons on a grid ...')
if ~isfield(options, 'dims')
    [axon_collection, dims] = setupAxonsGrid(original_axon_collection);
else
    [axon_collection, dims] = setupAxonsGrid(original_axon_collection, options.dims);
end
disp('done')

if ~isfield(options, 'mask')
    options.mask = zeros(dims);
    options.mask(round(dims(1)/3):round(2*dims(1)/3), round(dims(2)/3):round(2*dims(2)/3)) = 1;
end

% Axon packing
disp('process packing ...')
[axon_collection, FVF_packed_model] = packAxons(axon_collection, options.mask, options.max_iteration, options.max_FVF, options.packing_speed, options.plot_model);
disp('done')
disp(['FVF packed model : ' num2str(FVF_packed_model)]);

% Axon dispersion
if (options.dispersion_mode == 'remove')
    disp('remove axons ...')
    [axon_collection, FVF] = removeAxons(axon_collection, expected_FVF, options.tolerance, options.mask, options.plot_model);
elseif (options.dispersion_mode == 'spread')
    disp('spread axons ...')
    [axon_collection, FVF] = repulseAxons(axon_collection, expected_FVF, options.tolerance, options.mask, options.plot_model);
    axon_collection = avoidAxonOverlap(axon_collection, dims);
else 
    error('dispersion mode should be remove or spread');
end
disp(['current FVF : ' num2str(FVF)]);
disp('done')

axon_collection = convertAxonDataToRoundValues(axon_collection);
axon_collection_save = axon_collection;

% Change g-ratio (optional)
if options.expected_g_ratio
    disp('change g ratio ...')
    axon_collection = changeGRatio(axon_collection_save, options.expected_g_ratio, dims);
    disp('done')
end

[Model, ZoomedModel, FVF, g_ratio] = createModelFromData(axon_collection, options.mask, options.plot_model);

% Save model
mask = options.mask;
if options.save_model
    disp('save model ...')
    disp(options.save_model);
    save(options.save_model, 'Model', 'ZoomedModel', 'FVF', 'g_ratio', 'axon_collection', 'dims', 'mask')
    disp('done')
end

end