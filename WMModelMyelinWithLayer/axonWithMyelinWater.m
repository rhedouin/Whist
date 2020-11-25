function [simpleModel, model, mask] = axonWithMyelinWater(FVF, gRatio, nb_myelin_water_layer, myelin_water_proportion, N)

model = zeros(2*N + 1);

[xGrid, yGrid] = meshgrid(-N:N,-N:N);
dist = sqrt(xGrid.^2 + yGrid.^2);

mask = zeros(size(model));
mask(logical(dist < N)) = 1;

original_myelinated_axon_radius = N * sqrt(FVF);
intra_axonal_radius = gRatio * original_myelinated_axon_radius;
myelin_thickness = (1 - gRatio) * original_myelinated_axon_radius;

myelin_phospholipid_proportion = 1 - myelin_water_proportion;

if nb_myelin_water_layer ~= 0
    unit_myelin_phospholipid_layer_thickness = myelin_phospholipid_proportion * myelin_thickness / (nb_myelin_water_layer+1);
    unit_myelin_water_layer_thickness = myelin_water_proportion * myelin_thickness / nb_myelin_water_layer;
else
    unit_myelin_phospholipid_layer_thickness = 0;
    unit_myelin_water_layer_thickness = myelin_thickness;
end
    
unit_combine_layer_thickness = unit_myelin_phospholipid_layer_thickness + unit_myelin_water_layer_thickness;

model(dist < intra_axonal_radius) = 0.5; 
simpleModel = model;

simpleModel(logical((dist < intra_axonal_radius + myelin_thickness) .* (dist > intra_axonal_radius))) = 1; 
simpleModel(logical((dist > intra_axonal_radius + myelin_thickness) .* (dist < N))) = 0;
simpleModel(logical(dist > N)) = -1;

if nb_myelin_water_layer == 0
    model = simpleModel;
    return;
else
    
for k = 0:nb_myelin_water_layer-1
    model(logical((dist > intra_axonal_radius + k*unit_combine_layer_thickness) .* (dist < intra_axonal_radius + k*unit_combine_layer_thickness + unit_myelin_phospholipid_layer_thickness))) = 1; 
    model(logical((dist > intra_axonal_radius + k*unit_combine_layer_thickness + unit_myelin_phospholipid_layer_thickness) .* (dist < intra_axonal_radius + (k+1)*unit_combine_layer_thickness))) = 0.75; 
end

model(logical((dist > intra_axonal_radius + nb_myelin_water_layer*unit_combine_layer_thickness) .* (dist < intra_axonal_radius + nb_myelin_water_layer*unit_combine_layer_thickness + unit_myelin_phospholipid_layer_thickness))) = 1; 
model(logical((dist > intra_axonal_radius + myelin_thickness) .* (dist < N))) = 0;
model(logical(dist > N)) = -1;


end
