function [Model, ZoomedModel, FVF, g_ratio] = createModelFromData(axonCollection, mask, plot_model)

    dims = size(mask);
    [mask_x_index, mask_y_index] = find(mask);
    mask_min_x = min(mask_x_index);
    mask_max_x = max(mask_x_index);
    mask_min_y = min(mask_y_index);
    mask_max_y = max(mask_y_index);
    
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
        
    AVF = length(find(ZoomedModel == 0.5));
    MVF = length(find(ZoomedModel == 1));
    
    FVF = (AVF + MVF) / length(ZoomedModel(:));
    g_ratio = sqrt(AVF / (AVF + MVF));

    if plot_model
        imagesc(Model);
        hold on
        
        plotRectangle = 1;
        if plotRectangle
            rectangle('Position',[mask_min_x, mask_min_y, mask_max_x - mask_min_x, mask_max_y - mask_min_y],'EdgeColor', 'r', 'LineWidth', 3);
            title(['FVF = ' num2str(FVF) ', g-ratio = ' num2str(g_ratio)])
            set(gca, 'FontSize', 15)
        end
    end
end
