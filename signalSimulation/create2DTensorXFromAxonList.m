function [tensor_X, total_model, phimap] = create2DTensorXFromAxonList(axonlist,model_params)

total_X = zeros(model_params.dims(1),model_params.dims(2),3,3);

myelin_Xi = [model_params.myelin.xi 0 0; 0 model_params.myelin.xi 0; 0 0 model_params.myelin.xi];
myelin_Xa = [model_params.myelin.xa 0 0; 0 -model_params.myelin.xa/2 0; 0 0 -model_params.myelin.xa/2];

if isfield(model_params.intra_axonal, 'xi')
    intra_axonal_Xi = [model_params.intra_axonal.xi 0 0; 0 model_params.intra_axonal.xi 0; 0 0 model_params.intra_axonal.xi];
end

if isfield(model_params.extra_axonal, 'xi')
    extra_axonal_Xi = [model_params.extra_axonal.xi 0 0; 0 model_params.extra_axonal.xi 0; 0 0 model_params.extra_axonal.xi];
    total_X = permute(repmat(extra_axonal_Xi, [1 1 model_params.dims]), [3 4 1 2]);
end

total_model = zeros(model_params.dims);
phimap = zeros(model_params.dims);
compute_phimap = 1;

grad_threshold = 0.0001;

for j = 1:length(axonlist)

    map = zeros(model_params.dims);
    sub_myelin = axonlist(j).data;
    sub_intra_axonal = myelin2axon(sub_myelin);

    ind_myelin = sub2ind(model_params.dims,sub_myelin(:,1),sub_myelin(:,2));  
    ind_intra_axonal = sub2ind(model_params.dims,sub_intra_axonal(:,1),sub_intra_axonal(:,2));    
 
    min_sub_myelin = min(sub_myelin);
    max_sub_myelin = max(sub_myelin);
    
    extra_space = 10;   
    small_map_model_params.dims = max_sub_myelin - min_sub_myelin + 2*extra_space;
    
    new_sub_myelin = sub_myelin - min_sub_myelin + extra_space;
    new_sub_intra_axonal = myelin2axon(new_sub_myelin);

    new_ind_myelin = sub2ind(small_map_model_params.dims,new_sub_myelin(:,1),new_sub_myelin(:,2));  
    new_ind_axon = sub2ind(small_map_model_params.dims,new_sub_intra_axonal(:,1),new_sub_intra_axonal(:,2));    
    
    small_map = zeros(small_map_model_params.dims);
    small_map(new_ind_myelin) = 1;
    small_map(new_ind_axon) = 2;
    
    map(ind_myelin) = 1;
    map(ind_intra_axonal) = 2;   

    total_model = total_model + map;

    sigma=2;
    
    smooth_small_map = imgaussfilt(small_map,sigma, 'FilterSize',5);   
    [gradient_magnitude,gradient_direction] = imgradient(smooth_small_map);
    
    % Counting variable for smoothing process
    c=1;
    while sum(gradient_magnitude(new_ind_myelin) <= grad_threshold) > 0 && c<20
        smooth_small_map = imgaussfilt(smooth_small_map,sigma, 'FilterSize',5);
        [gradient_magnitude, gradient_direction] = imgradient(smooth_small_map);
        c=c+1;
    end
    
    gradient_direction = (pi/180)*(gradient_direction + 90);

    for k = 1:size(new_sub_myelin,1)
        phi = mod(gradient_direction(new_sub_myelin(k,1),new_sub_myelin(k,2)) + 2*pi, 2*pi) - pi;

        R = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
        total_X(sub_myelin(k,1),sub_myelin(k,2),:,:) = myelin_Xi + R*myelin_Xa*inv(R);
        
        if compute_phimap
            phimap(sub_myelin(k,1),sub_myelin(k,2)) = angle( exp(1i*phi) );
        end
    end
    
    if isfield(model_params.intra_axonal, 'xi')
        for k = 1:size(sub_intra_axonal,1)
            total_X(sub_intra_axonal(k,1),sub_intra_axonal(k,2),:,:) = intra_axonal_Xi;
        end
    end
end

total_model = min(total_model, 2);
total_model(find(total_model)) = 1./total_model(find(total_model));

tensor_X(:,:,1) = total_X(:,:,1,1);
tensor_X(:,:,2) = total_X(:,:,1,2);
tensor_X(:,:,3) = total_X(:,:,1,3);
tensor_X(:,:,4) = total_X(:,:,2,2);
tensor_X(:,:,5) = total_X(:,:,2,3);
tensor_X(:,:,6) = total_X(:,:,3,3);

end
