function out = createModelFromAxons(N, FVF_expected, FVF_max, tol, suffix)

project_folder = '/project/3015069.04';
WMModels_folder = [project_folder '/WM_Models/'];
histology_Folder = [project_folder '/histology/axon_dictionary/';
 
FVF_num = ['FVF' num2str(100*FVF_expected) '_N' num2str(N) suffix  '/'];

iter_max = 15*N;    

disp('load axon dictionary ...')
load([axonFolder 'axonLargeDict.mat']);
disp('done')

listAxons = randi(length(axonDico), N, 1);
shortAxonDico = axonDico(listAxons);

% packing process of the axons
side = 1000;
disp('setupAxonGrid ...')
[axon_collection, side] = setupAxonsGrid(shortAxonDico,side);
disp('done')

disp('process packing ...')
[pts, axon_shapes, TotalModel, ZoomedModel, FVF] = packAxons(x0, axon_shapes, side, iter_max, FVF_max);
disp('done')

disp('remove axons ...')
[TotalModel, ZoomedModel, pts, axon_shapes, currentFVF, listAxons] = removeAxons(pts, axon_shapes, FVF_expected, tol, side, shortAxonDico);
disp('done')


end





