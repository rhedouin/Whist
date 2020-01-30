function [axon_collection, FVF_current] = repulseAxons(axon_collection, FVF_expected, tol, dims)

pts = cat(1,axon_collection(:).Centroid);

plot = 1;
step = 0.01;

[~, ~, FVF_current] = createModelFromData(axon_collection, dims, plot);

N = length(axon_collection);

iter = 1;
while or((FVF_current < FVF_expected - tol),(FVF_current > FVF_expected))
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
    
    [~, ~, FVF_current] = createModelFromData(axon_collection, dims, plot);

    if (FVF_current < FVF_expected - tol)
        axon_collection = axon_collection_old;
        step = step/2;
    end
    
end

end