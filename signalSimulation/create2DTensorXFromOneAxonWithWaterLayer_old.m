function [tensor_X, total_model, phimap] = create2DTensorXFromOneAxonWithWaterLayer(model, simpleModel, model_parameters)
% Create 2D susceptibility tensor (see Tianyou Xu 2018)
model_parameters.dims = size(model);

total_X = zeros(model_parameters.dims(1),model_parameters.dims(2),3,3);

myelin_Xi = [model_parameters.myelin.xi 0 0; 0 model_parameters.myelin.xi 0; 0 0 model_parameters.myelin.xi];
myelin_Xa = [model_parameters.myelin.xa 0 0; 0 -model_parameters.myelin.xa/2 0; 0 0 -model_parameters.myelin.xa/2];

if (isfield(model_parameters, 'intra_axonal') && isfield(model_parameters.intra_axonal, 'xi'))
    intra_axonal_Xi = [model_parameters.intra_axonal.xi 0 0; 0 model_parameters.intra_axonal.xi 0; 0 0 model_parameters.intra_axonal.xi];
end

if (isfield(model_parameters, 'extra_axonal') && isfield(model_parameters.extra_axonal, 'xi'))
    extra_axonal_Xi = [model_parameters.extra_axonal.xi 0 0; 0 model_parameters.extra_axonal.xi 0; 0 0 model_parameters.extra_axonal.xi];
    total_X = permute(repmat(extra_axonal_Xi, [1 1 model_parameters.dims]), [3 4 1 2]);
end

if (isfield(model_parameters, 'myelin_water_layer') && isfield(model_parameters.myelin_water_layer, 'xi'))
    myelin_water_layer_Xi = [model_parameters.myelin_water_layer.xi 0 0; 0 model_parameters.myelin_water_layer.xi 0; 0 0 model_parameters.myelin_water_layer.xi];
end

total_model = zeros(model_parameters.dims);
phimap = zeros(model_parameters.dims);
compute_phimap = 1;

grad_threshold = 0.0001;

map = zeros(model_parameters.dims);

ind_intra_axon = find(model == 2);
ind_extra_axon = find(model == 1);
ind_myelin = find(model == 0);
ind_myelin_water = find(model == 0.5);

[sub_intra_axonal_x, sub_intra_axonal_y] = ind2sub(model_parameters.dims, ind_intra_axon);
[sub_extra_axonal_x, sub_extra_axonal_y] = ind2sub(model_parameters.dims, ind_extra_axon);
[sub_myelin_x, sub_myelin_y] = ind2sub(model_parameters.dims, ind_myelin);
[sub_myelin_water_x, sub_myelin_water_y] = ind2sub(model_parameters.dims, ind_myelin_water);

min_sub_myelin = [min(sub_myelin_x), min(sub_myelin_y)];
max_sub_myelin = [max(sub_myelin_x), max(sub_myelin_y)] ;

sigma=2;

smooth_map = imgaussfilt(simpleModel, sigma, 'FilterSize',5);
[gradient_magnitude,gradient_direction] = imgradient(smooth_map);

% Counting variable for smoothing process
c=1;
while sum(gradient_magnitude(ind_myelin) <= grad_threshold) > 0 && c<200
    smooth_map = imgaussfilt(smooth_map, sigma, 'FilterSize',7);
    [gradient_magnitude, gradient_direction] = imgradient(smooth_map);
    c = c+1
    if  (mod(c, 40) == 0)
        figure
        imagesc(smooth_map)
    end
end

if (isfield(model_parameters, 'intra_axonal') && isfield(model_parameters.intra_axonal, 'xi'))
    for k = 1:size(sub_intra_axonal_x,1)
        total_X(sub_intra_axonal_x(k), sub_intra_axonal_y(k),:,:) = intra_axonal_Xi;
    end   
end

if (isfield(model_parameters, 'extra_axonal') && isfield(model_parameters.extra_axonal, 'xi'))
    for k = 1:size(sub_extra_axonal_x,1)
        total_X(sub_extra_axonal_x(k), sub_extra_axonal_y(k),:,:) = extra_axonal_Xi;
    end   
end

if (isfield(model_parameters, 'myelin_water_layer') && isfield(model_parameters.myelin_water_layer, 'xi'))
    for k = 1:size(sub_myelin_water_x,1)
        total_X(sub_myelin_water_x(k), sub_myelin_water_y(k),:,:) = myelin_water_layer_Xi;
    end   
end

gradient_direction = (pi/180)*(gradient_direction + 90);

for k = 1:size(sub_myelin_x,1)
    phi = mod(gradient_direction(sub_myelin_x(k), sub_myelin_y(k)) + 2*pi, 2*pi) - pi;
    
    R = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    total_X(sub_myelin_x(k), sub_myelin_y(k),:,:) = myelin_Xi + R*myelin_Xa*inv(R);
    
    if compute_phimap
        phimap(sub_myelin_x(k), sub_myelin_y(k)) = angle(exp(1i*phi));
    end
end

figure;
imagesc(phimap)

keyboard;
total_model = min(total_model, 2);
total_model(find(total_model)) = 1./total_model(find(total_model));

tensor_X(:,:,1) = total_X(:,:,1,1);
tensor_X(:,:,2) = total_X(:,:,1,2);
tensor_X(:,:,3) = total_X(:,:,1,3);
tensor_X(:,:,4) = total_X(:,:,2,2);
tensor_X(:,:,5) = total_X(:,:,2,3);
tensor_X(:,:,6) = total_X(:,:,3,3);

if isfield(model_parameters, 'no_mask_tensor_map') && model_parameters.no_mask_tensor_map == 1
    mask_replic = repmat(model_parameters.mask,[1 1 6]);
    tensor_X(~mask_replic) = 0;
end

end
