function axon_collection = createOne2DWMModel(expected_FVF, expected_g_ratio, suffix, mode)

% parameters  
side = 1000;
dims = [side side];

iter_max = 5000;
maxFVF = 0.85;
plot = 1;

dico_folder = '/project/3015069.01/model/data/dicos/WM_Models/N400/';
 
FVF_num = ['FVF' num2str(100*expected_FVF) '_N400' suffix  '/'];
FVF_folder = [dico_folder FVF_num];
mkdir(FVF_folder);

if ~exist(mode)
    mode = 'spread';
end

axonFolder = '/project/3015069.04/histology/axon_dictionary/';
 
disp('load axon dictionary ...')
load([axonFolder 'axonLargeDict.mat']);
disp('done')

list_axons = randi(length(axonDico), 400, 1);
original_axon_collection = axonDico(list_axons);

for k = 1:length(original_axon_collection)
    original_axon_collection(k).data = double(original_axon_collection(k).data);
end

disp('setupAxonGrid ...')
[axon_collection, dims] = setupAxonsGrid(original_axon_collection,dims);
disp('done')

disp('process packing ...')
[axon_collection, FVF_packed_model] = packAxons(axon_collection, dims, iter_max, maxFVF);
disp('done')
disp(['FVF_packed_model : ' num2str(FVF_packed_model)]);

tol = 0.1;
if mode == 'remove'
    disp('remove axons ...')
    [axon_collection, FVF_current] = removeAxons(axon_collection, expected_FVF, tol, dims);
elseif mode == 'spread'
    disp('spread axons ...')
    [axon_collection, FVF_current] = repulseAxons(axon_collection, expected_FVF, tol, dims);
    axon_collection = avoidAxonOverlap(axon_collection, dims);
else 
    error('mode should be remove or spread');
end
disp('done')
disp(['current FVF : ' num2str(FVF_current)]);

[Model, ZoomedModel, current_FVF, ~] = createModelFromData(axon_collection, dims, plot);

axon_collection = convertAxonDataToRoundValues(axon_collection);

disp('save model ...')
save([dico_folder FVF_num 'AxonMap_FVF' num2str(100*expected_FVF) suffix '_N400.mat'], 'Model', 'ZoomedModel', 'current_FVF', 'axon_collection', 'dims')
disp('done')

axon_collection_save = axon_collection;

for g_ratio = expected_g_ratio
    
    disp('change g ratio ...')
    [axon_collection, current_g_ratio] = changeGRatio(axon_collection_save, g_ratio, dims);
    disp('done')
    
    [Model, ZoomedModel, current_FVF, ~] = createModelFromData(axon_collection, dims, plot);
    
    disp('save model ...')
    save([dico_folder FVF_num 'AxonMap_FVF' num2str(100*expected_FVF) suffix '_gratio_' num2str(100*g_ratio) '_N400.mat'], 'Model', 'ZoomedModel', 'current_FVF', 'current_g_ratio', 'axon_collection', 'dims')
    disp('done')
    
end
end