function [tensor_X, phimap] = create2DTensorXFromOneAxon(model, model_parameters)
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


phimap = zeros(model_parameters.dims);
compute_phimap = 1;

ind_intra_axon = find(model == 0.5);
ind_extra_axon = find(model == 0);
ind_myelin = find(model == 1);

[sub_intra_axonal_x, sub_intra_axonal_y] = ind2sub(model_parameters.dims, ind_intra_axon);
[sub_extra_axonal_x, sub_extra_axonal_y] = ind2sub(model_parameters.dims, ind_extra_axon);
[sub_myelin_x, sub_myelin_y] = ind2sub(model_parameters.dims, ind_myelin);

N = (model_parameters.dims(1)-1)/2;

[xGrid, yGrid] = meshgrid(-N:N,-N:N);
dist = sqrt(xGrid.^2 + yGrid.^2);

[~, gradient_direction] = imgradient(dist);

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

if (isfield(model_parameters, 'myelin') && isfield(model_parameters.myelin, 'xi'))
    for k = 1:size(sub_myelin_x,1)
        total_X(sub_myelin_x(k), sub_myelin_y(k),:,:) = myelin_Xi;
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

tensor_X(:,:,1) = total_X(:,:,1,1);
tensor_X(:,:,2) = total_X(:,:,1,2);
tensor_X(:,:,3) = total_X(:,:,1,3);
tensor_X(:,:,4) = total_X(:,:,2,2);
tensor_X(:,:,5) = total_X(:,:,2,3);
tensor_X(:,:,6) = total_X(:,:,3,3);

end
