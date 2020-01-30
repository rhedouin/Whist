function [Model, ZoomedModel, FVF, g_ratio] = createModelFromData(axonCollection, dims, plot)

    Model = zeros(dims);
    allMyelin = [];
    allAxon = [];

    for k= 1:length(axonCollection)
        myelin{k} = round(axonCollection(k).data);
        axon{k} = myelin2axon(myelin{k});
        allMyelin = [allMyelin; myelin{k}];
        allAxon = [allAxon; axon{k}];      
    end
    
    indexMyelin = sub2ind(dims, allMyelin(:,1), allMyelin(:,2)); 
    indexAxon = sub2ind(dims, allAxon(:,1), allAxon(:,2));

    Model(indexMyelin) = 1;
    Model(indexAxon)   = 0.5;
  
    ZoomedModel = Model(round(dims(1)/3):round(2*dims(1)/3), round(dims(2)/3):round(2*dims(2)/3));
        
    AVF = length(find(ZoomedModel == 0.5));
    MVF = length(find(ZoomedModel == 1));
    
    FVF = (AVF + MVF) / length(ZoomedModel(:));
    g_ratio = sqrt(AVF / (AVF + MVF));

    if plot
        imagesc(Model);
        hold on
        
        plotRectangle = 1;
        if plotRectangle
            rectangle('Position',[round(dims(1)/3), round(dims(1)/3), round(dims(2)/3), round(dims(2)/3)],'EdgeColor', 'r', 'LineWidth', 3);
            title(['FVF = ' num2str(FVF)])
            set(gca, 'FontSize', 15)
        end
    end
end
