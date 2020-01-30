function [tensor_X, total_model, phimap] = createTensorXFromAxonList(axonlist,dims,xa,xi)

total_X = zeros(dims(1),dims(2),3,3);

Xi = [xi 0 0; 0 xi 0; 0 0 xi];
Xa = [xa 0 0; 0 -xa/2 0; 0 0 -xa/2];

total_model = zeros(dims);
phimap = zeros(dims);
compute_phimap = 1;

grad_threshold = 0.0001;

for j = 1:length(axonlist)
    %% From Wharton 12        
    map = zeros(dims);
    ind_myelin = axonlist(j).data;
    ind_axon = myelin2axon(ind_myelin);

    sub_myelin = sub2ind(dims,ind_myelin(:,1),ind_myelin(:,2));  
    sub_axon = sub2ind(dims,ind_axon(:,1),ind_axon(:,2));    
 
    min_ind_myelin = min(ind_myelin);
    max_ind_myelin = max(ind_myelin);
    
    extra_space = 10;   
    small_map_dims = max_ind_myelin - min_ind_myelin + 2*extra_space;
    
    new_ind_myelin = ind_myelin - min_ind_myelin + extra_space;
    new_ind_axon = myelin2axon(new_ind_myelin);

    new_sub_myelin = sub2ind(small_map_dims,new_ind_myelin(:,1),new_ind_myelin(:,2));  
    new_sub_axon = sub2ind(small_map_dims,new_ind_axon(:,1),new_ind_axon(:,2));    
    
    small_map = zeros(small_map_dims);
    small_map(new_sub_myelin) = 1;
    small_map(new_sub_axon) = 2;
    
    map(sub_myelin) = 1;
    map(sub_axon) = 2;   

    total_model = total_model + map;


    sigma=2;
    
    smooth_small_map = imgaussfilt(small_map,sigma, 'FilterSize',5);   
    [gradient_magnitude,gradient_direction] = imgradient(smooth_small_map);
    c=1; 
    while sum(gradient_magnitude(new_sub_myelin) <= grad_threshold) > 0 && c<20
        smooth_small_map = imgaussfilt(smooth_small_map,sigma, 'FilterSize',5);  
        [gradient_magnitude, gradient_direction] = imgradient(smooth_small_map);
        c=c+1;
    end
    
    gradient_direction = (pi/180)*(gradient_direction + 90);

    for k = 1:size(new_ind_myelin,1)
        phi = mod(gradient_direction(new_ind_myelin(k,1),new_ind_myelin(k,2)) + 2*pi, 2*pi) - pi;

        R = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
        total_X(ind_myelin(k,1),ind_myelin(k,2),:,:) = Xi+R*Xa*inv(R);
        
        if compute_phimap
            phimap(ind_myelin(k,1),ind_myelin(k,2)) = angle( exp(1i*phi) );
        end
    end
end

% Temporary solution to avoid overlad between myelin and intra-axon, need to be improve
total_model = min(total_model, 2);
total_model(find(total_model)) = 1./total_model(find(total_model));

tensor_X(:,:,1) = total_X(:,:,1,1);
tensor_X(:,:,2) = total_X(:,:,1,2);
tensor_X(:,:,3) = total_X(:,:,1,3);
tensor_X(:,:,4) = total_X(:,:,2,2);
tensor_X(:,:,5) = total_X(:,:,2,3);
tensor_X(:,:,6) = total_X(:,:,3,3);

end
