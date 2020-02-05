function [axon_collection, FVF_current] = repulseAxons(axon_collection, FVF_expected, tol, mask, plot_model)

pts = cat(1,axon_collection(:).Centroid);
dims = size(mask);

step = 0.01;

[~, ~, FVF_current] = createModelFromData(axon_collection, mask, plot_model);
iter_total = 0;

while or((FVF_current < FVF_expected - tol),(FVF_current > FVF_expected))
    iter_total = iter_total + 1;

    display(['currentFVF : ' num2str(FVF_current)]);

    axon_collection_old = axon_collection;
    pts_directions = (pts - [dims(1)/2 dims(2)/2]) ;
    pts = pts + step*pts_directions;
    
    k = 1;
    while k < length(axon_collection)
        axon_collection(k).data = axon_collection(k).data - axon_collection(k).Centroid + pts(k,:);
        
        % Remove axons 
        if ((min(pts(k, :)) < 100) || (max(pts(k, :)) > 900))
            axon_collection(k) = [];  
            pts(k,:) = [];
        else
            axon_collection(k).Centroid = pts(k,:);
            k = k + 1;
        end
        
    end
    
    if (mod(iter_total, 5) == 0)
        [~, ~, FVF_current] = createModelFromData(axon_collection, mask, plot_model);
        pause(0.01);
    else
        [~, ~, FVF_current] = createModelFromData(axon_collection, mask, 0);
    end

    if (FVF_current < FVF_expected - tol)
        axon_collection = axon_collection_old;
        [~, ~, FVF_current] = createModelFromData(axon_collection, mask, 0);
        step = step/2;
    end  
end

end