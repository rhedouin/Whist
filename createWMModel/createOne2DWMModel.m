function axon_collection = createOne2DWMModel(axon_dictionary_path, expected_FVF, options)

% set default parameters  
if ~isfield(options, 'max_FVF')
    options.max_FVF = 0.85;
end

if ~isfield(options, 'max_iteration')
    options.max_iteration = 5000;
end

if ~isfield(options, 'packing_speed')
    options.packing_speed = 0.5;
end

if ~isfield(options, 'side')
    options.side = 1000;
end

if ~isfield(options, 'dispersion_mode')
    options.dispersion_mode = 'spread';
end
 
if ~isfield(options, 'number_of_axons')
    options.number_of_axons = 400;
end

if ~isfield(options, 'plot')
    options.plot = 1;
end

if ~isfield(options, 'tolerance')
    options.tolerance = 0.1;
end

if ~isfield(options, 'mask')
    options.mask = zeros(dims);
    options.mask(round(dims(1)/3):round(2*dims(1)/3), round(dims(2)/3):round(2*dims(2)/3));
end

disp('load axon dictionary ...')
load(axon_dictionary_path);
disp('done')

list_axons = randi(length(axonDico), options.number_of_axons, 1);
original_axon_collection = axonDico(list_axons);

for k = 1:length(original_axon_collection)
    original_axon_collection(k).data = double(original_axon_collection(k).data);
end

disp('setup axons grid ...')
dims = [options.side options.side];
[axon_collection, dims] = setupAxonsGrid(original_axon_collection, dims);
disp('done')

disp('process packing ...')
[axon_collection, FVF_packed_model] = packAxons(axon_collection, options.mask, options.max_iteration, options.max_FVF, options.packing_speed);
disp('done')
disp(['FVF packed model : ' num2str(FVF_packed_model)]);

if (options.dispersion_mode == 'remove')
    disp('remove axons ...')
    [axon_collection, FVF] = removeAxons(axon_collection, expected_FVF, options.tolerance, options.mask);
elseif (options.dispersion_mode == 'spread')
    disp('spread axons ...')
    [axon_collection, FVF] = repulseAxons(axon_collection, expected_FVF, options.tolerance, options.mask);
    axon_collection = avoidAxonOverlap(axon_collection, dims);
else 
    error('dispersion mode should be remove or spread');
end

disp('done')
disp(['current FVF : ' num2str(FVF)]);

axon_collection = convertAxonDataToRoundValues(axon_collection);
axon_collection_save = axon_collection;

if options.expected_g_ratio
    disp('change g ratio ...')
    axon_collection = changeGRatio(axon_collection_save, options.expected_g_ratio, dims);
    disp('done')
end

[Model, ZoomedModel, FVF, g_ratio] = createModelFromData(axon_collection, options.mask, options.plot);

mask = options.mask;
if options.save_model
    disp('save model ...')
    save(options.save_model, 'Model', 'ZoomedModel', 'FVF', 'g_ratio', 'axon_collection', 'dims', 'mask')
    disp('done')
end

end