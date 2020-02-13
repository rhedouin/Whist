function [Model, ZoomedModel, FVF, g_ratio] = createModelFromData(axonCollection, mask, plot_model)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the WM model from the axonCollection.
%
% %%%%%%%%%%%%%% Inputs
% %%%% axon_collection is a structure of the white matter model where each
% element represent an axon with 
% % - data: corresponds to the myelin sheath
% % - intraAxon:  corresponds to the intra axonal (only required for 3D)
%
% %%%% The mask need to correspond to the axonCollection as it also
% contains the dimension of the model
% If you plot the model, it is represented with the red rectangle for 2D models and with a
% different color for 3D models (as it is more common to have a complex mask
% shape for a 3D model to ensure a FVF high enough)
%
% %%%%%%%% Outputs
% ZoomedModel represents the model within the mask area, it is only relevant if
% the mask is a rectangle or a cuboid
% The fiber volume fraction (FVF) and g-ratio are computed within the mask

%%%%%%%%%%%%%%%%%%%%%%%%%%%

dims = size(mask);
number_dims = ndims(mask);

if number_dims == 2
    [mask_x_sub, mask_y_sub] = find(mask);
    mask_min_x = min(mask_x_sub);
    mask_max_x = max(mask_x_sub);
    mask_min_y = min(mask_y_sub);
    mask_max_y = max(mask_y_sub);
    
    Model = zeros(dims);
    allMyelin = [];
    allAxon = [];
    
    for k= 1:length(axonCollection)
        myelin{k} = round(axonCollection(k).data);
        axon{k} = myelin2axon(myelin{k});
        allMyelin = [allMyelin; myelin{k}];
        allAxon = [allAxon; axon{k}];
    end
    
    myelin_index = sub2ind(dims, allMyelin(:,1), allMyelin(:,2));
    axon_index = sub2ind(dims, allAxon(:,1), allAxon(:,2));
    
    Model(myelin_index) = 1;
    Model(axon_index)   = 0.5;
    
    ZoomedModel = Model(mask_min_x:mask_max_x, mask_min_y:mask_max_y);
    
    AVF = length(find(Model.*mask == 0.5));
    MVF = length(find(Model.*mask == 1));
    
    FVF = (AVF + MVF) / sum(mask, 'all');
    g_ratio = sqrt(AVF / (AVF + MVF));
    if plot_model
        imagesc(Model);
        hold on
        
        plotRectangle = 1;
        if plotRectangle
            rectangle('Position',[mask_min_x, mask_min_y, mask_max_x - mask_min_x, mask_max_y - mask_min_y],'EdgeColor', 'r', 'LineWidth', 3);
            title(['WM model, FVF = ' num2str(FVF) ', g-ratio = ' num2str(g_ratio)])
            set(gca, 'FontSize', 12)
        end
    end

elseif number_dims == 3
    
    mask_ind = find(mask);
    
    [mask_x_sub, mask_y_sub, mask_z_sub] = ind2sub(dims, mask_ind);
    mask_min_x = min(mask_x_sub);
    mask_max_x = max(mask_x_sub);
    mask_min_y = min(mask_y_sub);
    mask_max_y = max(mask_y_sub);
    mask_min_z = min(mask_z_sub);
    mask_max_z = max(mask_z_sub);
    
    Model = zeros(dims);
    allMyelin = [];
    allAxon = [];
    
    for k= 1:length(axonCollection)
        myelin{k} = round(axonCollection(k).data);
        axon{k} = round(axonCollection(k).intraAxon);
        allMyelin = [allMyelin; myelin{k}];
        allAxon = [allAxon; axon{k}];
    end
    
    myelin_index = sub2ind(dims, allMyelin(:,1), allMyelin(:,2), allMyelin(:,3));
    axon_index = sub2ind(dims, allAxon(:,1), allAxon(:,2), allAxon(:,3));
    
    Model(myelin_index) = 1;
    Model(axon_index)   = 0.5;
    Model(mask == 0) = -1;
    
    ZoomedModel = Model(mask_min_x:mask_max_x, mask_min_y:mask_max_y, mask_min_z:mask_max_z);
    
    AVF = length(find(Model.*mask == 0.5));
    MVF = length(find(Model.*mask == 1));
    
    FVF = (AVF + MVF) / sum(mask, 'all');
    g_ratio = sqrt(AVF / (AVF + MVF));
    
    if plot_model
        h = figure('Name', '3D WM model');
        position = [10 10 990 890];
        h.Position = position;           
        
        subplot(221)
        text(0, 0.5, {'3D WM model', 'slice view of the middle of each dimension', ...
            '3 compartments : intra axonal, myelin, extra axonal', 'the rest being outside of the mask', ...
            ['FVF = ' num2str(FVF) ], ['g-ratio = ' num2str(g_ratio)]});
        axis off
        
        subplot(222)
        imagesc(squeeze(Model(round(dims(1)/2), :,:)));
        xlabel('y')
        ylabel('z')
        title(['x slice number ' num2str(round(dims(1)/2))])
          set(gca, 'FontSize', 12)

        subplot(223)       
        imagesc(squeeze(Model(:, round(dims(2)/2), :)));
        xlabel('x')
        ylabel('z')
        title(['y slice number ' num2str(round(dims(2)/2))])
          set(gca, 'FontSize', 12)

        subplot(224)   
        imagesc(squeeze(Model(:, :, round(dims(3)/2))));
        xlabel('x')
        ylabel('y')
        title(['z slice number ' num2str(round(dims(3)/2))])  
        set(gca, 'FontSize', 12)

    end
    Model(mask == 0) = 0;
end
