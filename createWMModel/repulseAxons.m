function [axon_collection, FVF_current] = repulseAxons(axon_collection, FVF_expected, tol, mask, plot_model)

[~, ~, FVF_current] = createModelFromData(axon_collection, mask, plot_model);

if FVF_current < FVF_expected
    error('The current FVF before axon dispersions is lower than the expected FVF and cannot be reach')
end

% Define model boundaries 
maxSize = 0;

N = length(axon_collection);
for k = 1:N
    currentSize = max(max(axon_collection(k).data)  - min(axon_collection(k).data ));
    if (currentSize > maxSize)
        maxSize = max(currentSize);
    end
end

dims = size(mask);
[mask_x_index, mask_y_index] = find(mask);
mask_min_x = min(mask_x_index);
mask_max_x = max(mask_x_index);
mask_min_y = min(mask_y_index);
mask_max_y = max(mask_y_index);

edge = 2;
remove_edge_min_x = max(mask_min_x - (maxSize + edge), (maxSize + edge));
remove_edge_max_x = min(mask_max_x + (maxSize + edge), dims(1) - (maxSize + edge));
remove_edge_min_y = max(mask_min_y - (maxSize + edge), (maxSize + edge));
remove_edge_max_y = min(mask_max_y + (maxSize + edge), dims(2) - (maxSize + edge));

pts = cat(1,axon_collection(:).Centroid);
dims = size(mask);

step = 0.01;

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
        
        % Remove axons out of boundaries    
        if ((min(pts(k, 1)) < remove_edge_min_x) || (max(pts(k, 1)) > remove_edge_max_x) || ...
                (min(pts(k, 2)) < remove_edge_min_y) || (max(pts(k, 2)) > remove_edge_max_y))
            
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