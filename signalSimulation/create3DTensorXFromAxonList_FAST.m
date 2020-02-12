function [TensorX, TotalModel,M_old] = createTensor3DXFromAxonList_FAST(axonlist,dims,xa,xi)


% Magnetic susceptibility of myelin (isotropic and anisotropic)
Xi = [xi 0 0; 0 xi 0; 0 0 xi];
Xa = [xa 0 0; 0 -xa/2 0; 0 0 -xa/2];

% Total_Model = zeros(dims);

% Preallocation
TotalModel = zeros(dims,'single');
TotalX = zeros(dims(1),dims(2),dims(3),3,3,'single');

phimap_3D = zeros(dims,'single');
thetamap_3D = zeros(dims,'single');

for j = 1:length(axonlist) 
    %% From Wharton 12
    
    % Counting variable for smoothing process
    b=0;

    map = zeros(dims,'single');
    ind_myelin = single(axonlist(j).data);
    ind_iA = single(axonlist(j).intraAxon);

    sub_myelin = sub2ind(dims,ind_myelin(:,1),ind_myelin(:,2),ind_myelin(:,3));
    if length(ind_iA(:)) <=3
        ind_iA(2,:) = ind_iA;
    else
    end
    
    sub_iA = sub2ind(dims,ind_iA(:,1),ind_iA(:,2),ind_iA(:,3)); 

    map(sub_myelin) = 1;
    map(sub_iA) = 2;
    
    min_ind_myelin = min(ind_myelin);
    max_ind_myelin = max(ind_myelin);
    min_ind_iA = min(ind_iA);
    
    extra_space = 10;   
    small_map_dims = max_ind_myelin - min_ind_myelin + 2*extra_space;
    
    new_ind_myelin = ind_myelin - min_ind_myelin + extra_space;
    new_ind_iA = ind_iA - min_ind_myelin + extra_space;

    new_sub_myelin = sub2ind(small_map_dims,new_ind_myelin(:,1),new_ind_myelin(:,2),new_ind_myelin(:,3)); 
    new_sub_iA = sub2ind(small_map_dims,new_ind_iA(:,1),new_ind_iA(:,2),new_ind_iA(:,3)); 
    
    small_map = zeros(small_map_dims,'single');
    small_map(new_sub_myelin) = 1;  

    M_old = regionprops3_1(small_map,'MajorAxis');

    small_map(new_sub_iA) = 2;

    clear new_sub_iA
    
    TotalModel = map+TotalModel;  
    clear map
    
    sigma= 2;

    smoothMap = imgaussfilt3(small_map, sigma, 'FilterSize', 5, 'Padding', 'replicate','FilterDomain', 'spatial'); %, 'Padding', 'replicate','FilterDomain', 'spatial'

    % Computes gradient magintude, phi and elevation 
    [Gmag, Gphi, Gelev] = imgradient3(smoothMap);

    while sum(round(Gmag(new_sub_myelin),1) <= 0.1) >0  && b<10
        smoothMap = imgaussfilt3(smoothMap, sigma, 'FilterSize', 5, 'Padding', 'replicate','FilterDomain', 'spatial'); %, 'Padding', 'replicate','FilterDomain', 'spatial'
        [Gmag, Gphi, Gelev] = imgradient3(smoothMap);
        b=b+1;
    end
    
    clear Gmag smoothMap b new_sub_myelin small_map
   
    phi = (pi/180)*(Gphi - 90);     
    clear Gphi
    Gelev_degree = (pi/180)*Gelev;     
    clear Gelev
    theta = Gelev_degree; 

    clear Gelev_degree

    for k = 1:size(new_ind_myelin,1)

        phi_rot = phi(new_ind_myelin(k,1), new_ind_myelin(k,2),new_ind_myelin(k,3));
        theta_rot = theta(new_ind_myelin(k,1), new_ind_myelin(k,2),new_ind_myelin(k,3));
        Rz = [cos(phi_rot) -sin(phi_rot) 0; sin(phi_rot) cos(phi_rot) 0; 0 0 1];
        Ry = [cos(theta_rot) 0 sin(theta_rot);0 1 0; -sin(theta_rot) 0 cos(theta_rot)];
        
        R = Rz*Ry;
        
        % Computation of total magnetic susceptibility X
        TotalX(ind_myelin(k,1),ind_myelin(k,2),ind_myelin(k,3),:,:) = Xi + R*Xa*inv(R);
        
        phimap_3D(ind_myelin(k,1),ind_myelin(k,2),ind_myelin(k,3)) = angle( exp(1i*phi_rot) );
        thetamap_3D(ind_myelin(k,1),ind_myelin(k,2),ind_myelin(k,3)) = angle( exp(1i*theta_rot) );
      
    end
clear phi_rot theta_rot R Rz Ry theta phi new_ind_myelin ind_myelin sub_iA sub_myelin

end
% keyboard;
clear map ind_iA axonlist new_ind_iA

% Creates values in model .5 & 1
TotalModel = min(TotalModel, 2);
TotalModel(find(TotalModel)) = 1./TotalModel(find(TotalModel));

% Magnetic susceptibility tensor
X1 = TotalX(:,:,:,1,1);
X2 = TotalX(:,:,:,1,2);
X3 = TotalX(:,:,:,1,3);
X4 = TotalX(:,:,:,2,2);
X5 = TotalX(:,:,:,2,3);
X6 = TotalX(:,:,:,3,3);

clear TotalX
TensorX = zeros([dims 6],'single');
% load('X1') use load function for LARGE models to avoid MATLAB to crash!
TensorX(:,:,:,1) = X1;
% clear X1
% load('X2')
TensorX(:,:,:,2) = X2;
% clear X2
% load('X3')
TensorX(:,:,:,3) = X3;
% clear X3
% load('X4')
TensorX(:,:,:,4) = X4;
% clear X4
% load('X5')
TensorX(:,:,:,5) = X5;
% clear X5
% load('X6')
TensorX(:,:,:,6) = X6;
% clear X6

end
