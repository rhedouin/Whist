function [tensor_X, total_model, axonlist] = create3DTensorXFromAxonList(axonlist, model_params)

total_X = zeros(model_params.dims(1),model_params.dims(2),model_params.dims(3),3,3,'single');

myelin_Xi = [model_params.myelin.xi 0 0; 0 model_params.myelin.xi 0; 0 0 model_params.myelin.xi];
myelin_Xa = [model_params.myelin.xa 0 0; 0 -model_params.myelin.xa/2 0; 0 0 -model_params.myelin.xa/2];

if isfield(model_params, 'intra_axonal') && isfield(model_params.intra_axonal, 'xi')
    intra_axonal_Xi = [model_params.intra_axonal.xi 0 0; 0 model_params.intra_axonal.xi 0; 0 0 model_params.intra_axonal.xi];
end

if isfield(model_params, 'extra_axonal') && isfield(model_params.extra_axonal, 'xi')
    extra_axonal_Xi = [model_params.extra_axonal.xi 0 0; 0 model_params.extra_axonal.xi 0; 0 0 model_params.extra_axonal.xi];
    total_X = permute(repmat(extra_axonal_Xi, [1 1 model_params.dims]), [3 4 5 1 2]);
end

% Preallocation
total_model = zeros(model_params.dims,'single');
 
% phimap_3D = zeros(model_params.dims,'single');
% thetamap_3D = zeros(model_params.dims,'single');

grad_threshold = 0.0001;

for j = 1:length(axonlist) 
    j
    %% From Wharton 12
    map = zeros(model_params.dims,'single');
    sub_myelin = single(axonlist(j).data);
    sub_intra_axonal = single(axonlist(j).intraAxon);

    ind_myelin = sub2ind(model_params.dims,sub_myelin(:,1),sub_myelin(:,2),sub_myelin(:,3));
    if length(sub_intra_axonal(:)) <=3
        sub_intra_axonal(2,:) = sub_intra_axonal;
    else
    end
    
    ind_intra_axonal = sub2ind(model_params.dims,sub_intra_axonal(:,1),sub_intra_axonal(:,2),sub_intra_axonal(:,3)); 

    map(ind_myelin) = 1;
    map(ind_intra_axonal) = 2;
    
    min_sub_myelin = min(sub_myelin);
    max_sub_myelin = max(sub_myelin);
    
    extra_space = 10;   
    small_map_dims = max_sub_myelin - min_sub_myelin + 2*extra_space;
    
    new_sub_myelin = sub_myelin - min_sub_myelin + extra_space;
    new_sub_intra_axonal = sub_intra_axonal - min_sub_myelin + extra_space;

    new_ind_myelin = sub2ind(small_map_dims,new_sub_myelin(:,1),new_sub_myelin(:,2),new_sub_myelin(:,3)); 
    new_ind_intra_axonal = sub2ind(small_map_dims,new_sub_intra_axonal(:,1),new_sub_intra_axonal(:,2),new_sub_intra_axonal(:,3)); 
    
    small_map = zeros(small_map_dims,'single');
    small_map(new_ind_myelin) = 1;  

    M_old = regionprops3_1(small_map,'MajorAxis');
    
    axonlist(j).majorAxis = M_old.MajorAxis;
    
    small_map(new_ind_intra_axonal) = 2;

    clear new_sub_intra_axonal
    
    total_model = map+total_model;  
    clear map
    
    sigma= 2;

    smooth_small_map = imgaussfilt3(small_map, sigma, 'FilterSize', 5, 'Padding', 'replicate','FilterDomain', 'spatial'); 
    % Computes gradient magintude, phi and elevation 
    [gradient_magnitude, gradient_phi, gradient_elevation] = imgradient3(smooth_small_map);

    % Counting variable for smoothing process
    c=1;  
    while sum(gradient_magnitude(new_ind_myelin) <= grad_threshold) > 0 && c<20
        smooth_small_map = imgaussfilt3(smooth_small_map, sigma, 'FilterSize', 5, 'Padding', 'replicate','FilterDomain', 'spatial'); 
        [gradient_magnitude, gradient_phi, gradient_elevation] = imgradient3(smooth_small_map);
        c=c+1;
    end
    
    clear gradient_magnitude smoothMap b small_map
   
    phi = (pi/180)*(gradient_phi - 90);     
    clear gradient_phi
    gradient_elevation_degree = (pi/180)*gradient_elevation;     
    clear gradient_elevation
    theta = gradient_elevation_degree; 

    clear Gelev_degree
    
    for k = 1:size(new_sub_myelin,1)
        
        phi_rot = phi(new_sub_myelin(k,1), new_sub_myelin(k,2),new_sub_myelin(k,3));
        theta_rot = theta(new_sub_myelin(k,1), new_sub_myelin(k,2),new_sub_myelin(k,3));
        Rz = [cos(phi_rot) -sin(phi_rot) 0; sin(phi_rot) cos(phi_rot) 0; 0 0 1];
        Ry = [cos(theta_rot) 0 sin(theta_rot);0 1 0; -sin(theta_rot) 0 cos(theta_rot)];
        
        R = Rz*Ry;
        
        % Computation of total magnetic susceptibility X
        total_X(sub_myelin(k,1),sub_myelin(k,2),sub_myelin(k,3),:,:) = myelin_Xi + R*myelin_Xa*inv(R);
        
%         phimap_3D(sub_myelin(k,1),sub_myelin(k,2),sub_myelin(k,3)) = angle( exp(1i*phi_rot) );
%         thetamap_3D(sub_myelin(k,1),sub_myelin(k,2),sub_myelin(k,3)) = angle( exp(1i*theta_rot) );     
    end
    
    if isfield(model_params, 'intra_axonal') && isfield(model_params.intra_axonal, 'xi')
        for k = 1:size(sub_intra_axonal,1)
            total_X(sub_intra_axonal(k,1),sub_intra_axonal(k,2),sub_intra_axonal(k,3), :,:) = intra_axonal_Xi;
        end
    end
clear phi_rot theta_rot R Rz Ry theta phi new_ind_myelin ind_myelin sub_iA sub_myelin

end
% keyboard;
clear map ind_iA new_ind_iA

% Creates values in model .5 & 1
total_model = min(total_model, 2);
total_model(find(total_model)) = 1./total_model(find(total_model));

% Magnetic susceptibility tensor
X1 = total_X(:,:,:,1,1);
X2 = total_X(:,:,:,1,2);
X3 = total_X(:,:,:,1,3);
X4 = total_X(:,:,:,2,2);
X5 = total_X(:,:,:,2,3);
X6 = total_X(:,:,:,3,3);

clear TotalX
tensor_X = zeros([model_params.dims 6],'single');
% load('X1') use load function for LARGE models to avoid MATLAB to crash!
tensor_X(:,:,:,1) = X1;
% clear X1
% load('X2')
tensor_X(:,:,:,2) = X2;
% clear X2
% load('X3')
tensor_X(:,:,:,3) = X3;
% clear X3
% load('X4')
tensor_X(:,:,:,4) = X4;
% clear X4
% load('X5')
tensor_X(:,:,:,5) = X5;
% clear X5
% load('X6')
tensor_X(:,:,:,6) = X6;
% clear X6

end
