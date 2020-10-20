% test Axon Function 
close all;
clear;

FVF = 0.5;
gRatio = 0.7;
nbWaterLayer = 0;
myelinProportion = 0.8;
N = 50;

[simpleModel, model] = axonWithMyelinLayer(FVF, gRatio, nbWaterLayer, myelinProportion, N);
dims = size(simpleModel);

model_parameters.myelin.xi = -0.1;
model_parameters.myelin.xa = -0.1;

model_parameters.intra_axonal.xi = 0;
model_parameters.extra_axonal.xi = 0;
model_parameters.myelin_water_layer.xi = 0;

tensor_X = create2DTensorXFromOneAxonWithWaterLayer(model, simpleModel, model_parameters);

oneAxon.Centroid = [N, N];
listMyelin = find(simpleModel < 0.6);

[subMyelin_x, subMyelin_y] = ind2sub(dims, listMyelin);
oneAxon.data = [subMyelin_x, subMyelin_y];

model_parameters.dims = dims;

[tensor_X_simple, ~, ~] = create2DTensorXFromAxonList(oneAxon, model_parameters);

figure
subplot(231)
imagesc(tensor_X(:,:,1));
colorbar;

subplot(232)
imagesc(tensor_X(:,:,2));
colorbar;

subplot(233)
imagesc(tensor_X(:,:,3));
colorbar;

subplot(234)
imagesc(tensor_X(:,:,4));
colorbar;

subplot(235)
imagesc(tensor_X(:,:,5));
colorbar;

subplot(236)
imagesc(tensor_X(:,:,6));
colorbar;


figure
subplot(231)
imagesc(tensor_X_simple(:,:,1));
colorbar;

subplot(232)
imagesc(tensor_X_simple(:,:,2));
colorbar;

subplot(233)
imagesc(tensor_X_simple(:,:,3));
colorbar;

subplot(234)
imagesc(tensor_X_simple(:,:,4));
colorbar;

subplot(235)
imagesc(tensor_X_simple(:,:,5));
colorbar;

subplot(236)
imagesc(tensor_X_simple(:,:,6));
colorbar;