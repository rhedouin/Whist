function [axon_collection, model, zoomed_model] = createOne2DWMModel(axon_dictionary_path, model_params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Input required
%  axon_dictionary_path   : input axon dictionary path (provided in the
%  data folder)
 
% %%%%%%%%%% Input optional 
% model_params is a structure which contain all the algorithm options describe below
% The given values are the one set by defauts
%
% %%%% White matter model 
% model_params.number_of_axons = 400;
% % number of axons put on the original grid
% model_params.dims = estimated from the number of axons as well as the axon size 
% % dims of the original axon grid 
% model_params.mask = estimated from the original fiber volume fraction (FVF)
% % Represents the final WM model area after axon packing where the actual
% FVF, g-ratio are estimated
% 
% %%%% Axons packing 
% model_params.max_FVF = 0.85;
% % FVF value where the axon packing is stopped (it cannot be much higher than
% % 0.85 )
% model_params.max_iteration = 5000;
% % Max number of iteration of the axon packing (stop the axon packing if it
% % cannot reach the max_FVF)
% model_params.packing_speed = 0.5;
% % The packing speed weights the attraction/repulsion of the axons. A higher
% % value accelerates the packing process but can create axons overlap
% 
% %%%% Axons dispersion
% model_params.expected_FVF = 0.7; 
% %  FVF of the final WM model after remove or spread the axons
% model_params.dispersion_mode = 'spread'; 
% % Dispersion mode can be 
% % - remove where the axons are randomly remove from the packed model
% % - spread where the axons are spread from the grid center
% model_params.tolerance = 0.001;
% % Tolerance between expected FVF and actual FVF of the model
% 
% %%%% Change g-ratio
% model_params.expected_g_ratio = 0.6;
% % Change the myelin thickness to reach an expected g-ratio
% 
% %%%% Plot / save model
% model_params.plot_model = 1;
% model_params.save_model = '/project/3015069.04/WM_Models/toto.mat';
%
% %%%%%%%%%% Outputs
% %%%% axon_collection is a structure of the white matter model where each
% element represent an axon with one required field
% % - data which corresponds to the myelin sheath
% % - Centroid (required to run this function but not to simulate field
% perturbation)
% optional fields: gRatio, axonEquiDiameter, myelinThickness
%
% %%%% Model is the final WM model where 0 is extra axonal, 1 myelin, 0.5
% intra axonal
% %%%% Zoomed model is the model within the mask

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default parameters  
if ~exist('model_params')
    model_params.null = 0;
end

if ~isfield(model_params, 'max_FVF')
    model_params.max_FVF = 0.85;
end

if ~isfield(model_params, 'max_iteration')
    model_params.max_iteration = 5000;
end

if ~isfield(model_params, 'packing_speed')
    model_params.packing_speed = 0.5;
end

if ~isfield(model_params, 'dispersion_mode')
    model_params.dispersion_mode = 'spread';
end
 
if ~isfield(model_params, 'number_of_axons')
    model_params.number_of_axons = 400;
end

if ~isfield(model_params, 'plot_model')
    model_params.plot_model = 1;
end

if ~isfield(model_params, 'tolerance')
    model_params.tolerance = 0.1;
end

if  isfield(model_params, 'save_model')
    % makes sure directory exists
    a=max([find(model_params.save_model=='/'), find(model_params.save_model=='\')]);
    if ~isfolder(model_params.save_model(1:a))
        mkdir(model_params.save_model(1:a));
    end    
end

% Load dictionary
disp('load axon dictionary ...')
load(axon_dictionary_path);
disp('done')

disp(['randomly select ' num2str(model_params.number_of_axons) ' axon shapes ...']);
list_axons = randi(length(axonDico), model_params.number_of_axons, 1);
original_axon_collection = axonDico(list_axons);
for k = 1:length(original_axon_collection)
    original_axon_collection(k).data = double(original_axon_collection(k).data);
end
disp('done')

% Setup axons grid
disp('setup axons on a grid ...')
if ~isfield(model_params, 'dims')
    [axon_collection, dims] = setupAxonsGrid(original_axon_collection);
else
    [axon_collection, dims] = setupAxonsGrid(original_axon_collection, model_params.dims);
end
disp('done')

if ~isfield(model_params, 'mask')
    model_params.mask = createAdaptedMask(axon_collection, dims);
end

[model, zoomed_model, FVF, g_ratio] = createModelFromData(axon_collection, model_params.mask, model_params.plot_model);

% Axon packing
disp('process packing ...')
[axon_collection, FVF_packed_model] = packAxons(axon_collection, model_params.mask, model_params.max_iteration, model_params.max_FVF, model_params.packing_speed, model_params.plot_model);
disp('done')
disp(['FVF packed model : ' num2str(FVF_packed_model)]);

% Axon dispersion (optional)
if model_params.expected_FVF
    if (model_params.dispersion_mode == 'remove')
        disp('remove axons ...')
        [axon_collection, FVF] = removeAxons(axon_collection, model_params.expected_FVF, model_params.tolerance, model_params.mask, model_params.plot_model);
    elseif (model_params.dispersion_mode == 'spread')
        disp('spread axons ...')
        [axon_collection, FVF] = repulseAxons(axon_collection, model_params.expected_FVF, model_params.tolerance, model_params.mask, model_params.plot_model);
        axon_collection = avoidAxonOverlap(axon_collection, dims);
    else
        error('dispersion mode should be remove or spread');
    end
    disp(['current FVF : ' num2str(FVF)]);
    disp('done')
end

axon_collection = convertAxonDataToRoundValues(axon_collection);
axon_collection_save = axon_collection;

% Change g-ratio (optional)
if model_params.expected_g_ratio
    disp('change g ratio ...')
    axon_collection = changeGRatio(axon_collection_save, model_params.expected_g_ratio, model_params.mask);
    disp('done')
end

[model, zoomed_model, FVF, g_ratio] = createModelFromData(axon_collection, model_params.mask, model_params.plot_model);

% Save model
mask = model_params.mask;
if model_params.save_model
    disp('save model ...')
    disp(model_params.save_model);
    try
    save(model_params.save_model, 'model', 'zoomed_model', 'FVF', 'g_ratio', 'axon_collection', 'dims', 'mask')
    catch
        display ('failed to save...')
    end    
    disp('done')
end

end