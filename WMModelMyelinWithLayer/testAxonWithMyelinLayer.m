clear;
close all;

N = 200;
model = zeros(2*N + 1);

[xGrid, yGrid] = meshgrid(-N:N,-N:N);

dist = sqrt(xGrid.^2 + yGrid.^2);

imagesc(dist)

totalMyelinThickness = 150;
axonThickness = 30;
ratioMyelinWaterLayerThickness = 1;
nbWaterLayer = 20;
totalNbLayer = 2*nbWaterLayer + 1;
unitMyelinLayerThickness = (ratioMyelinWaterLayerThickness  / (ratioMyelinWaterLayerThickness + 1))* totalMyelinThickness / (totalNbLayer/2);
unitWaterLayerThickness = (1 / (ratioMyelinWaterLayerThickness + 1))* totalMyelinThickness / (totalNbLayer/2);
unitCombineLayerThickness = unitMyelinLayerThickness + unitWaterLayerThickness;

model(dist < axonThickness) = 2; 
simpleModel = model;

simpleModel(logical((dist < axonThickness + totalMyelinThickness) .* (dist > axonThickness))) = 0; 
simpleModel(logical(dist > axonThickness + totalMyelinThickness)) = 1;

figure
imagesc(simpleModel)

for k = 0:nbWaterLayer-1
    k
    model(logical((dist > axonThickness + k*unitCombineLayerThickness) .* (dist < axonThickness + k*unitCombineLayerThickness + unitMyelinLayerThickness))) = 0; 
    model(logical((dist > axonThickness + k*unitCombineLayerThickness + unitMyelinLayerThickness) .* (dist < axonThickness + (k+1)*unitCombineLayerThickness))) = 0.5; 
end

model(logical((dist > axonThickness + nbWaterLayer*unitCombineLayerThickness) .* (dist < axonThickness + nbWaterLayer*unitCombineLayerThickness + unitMyelinLayerThickness))) = 0; 
model((logical(dist > axonThickness + totalMyelinThickness))) = 1;

%model(logical((dist > axonThickness + totalMyelinThickness - unitMyelinLayerThickness) .* (dist < axonThickness + totalMyelinThickness))) = 5; 

figure
imagesc(model)
