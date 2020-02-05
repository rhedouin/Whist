function [axon_collection, mean_g_ratio] = changeGRatio(axon_collection, expected_g_ratio, dims)

g_list = cat(1,axon_collection(:).gRatio);
mean_g_ratio = mean(g_list);

pdf_list = cat(1,(axon_collection(:).axonEquivDiameter)).^2;

if (mean_g_ratio > expected_g_ratio)
    mode = 'extend';
elseif (mean_g_ratio < expected_g_ratio)
    mode = 'shrink';
else
    error('Expected GRatio');
end

stop = 0;
it = 0;

while stop == 0
    
    it = it+1;
    
    cdf_list = cumsum([0; pdf_list]);
    if (cdf_list(end) == 0)
        error('the expected gRatio could no be achieved');
    else
        cdf_list = cdf_list / cdf_list(end);
    end
    
    k = sum(rand >= cdf_list);
    
    pdf_list(k) = pdf_list(k)/2;
    
    clear current_myelin current_axon

    myelin_map = zeros(dims);
    
    old_myelin = round(axon_collection(k).data);    
    old_myelin_index = sub2ind(dims,old_myelin(:,1),old_myelin(:,2));
    myelin_map(old_myelin_index) = 1;
    
    old_axon = myelin2axon(old_myelin);
    old_axon_index = sub2ind(dims,old_axon(:,1),old_axon(:,2));
    
    nb_pixel = length([old_myelin_index; old_axon_index]);
    
    switch mode
        case 'extend'
            
            dilate_myelin_map = imdilate(myelin_map, ones(3,3));
            current_myelin_index = find(dilate_myelin_map);
            current_myelin_index = intersect(current_myelin_index,[old_myelin_index; old_axon_index]);
            
            [current_myelin(:,1), current_myelin(:,2)] = ind2sub(dims,current_myelin_index);
            
            current_axon = myelin2axon(current_myelin);
            current_axon_index = sub2ind(dims,current_axon(:,1),current_axon(:,2));
            
            if (length(current_myelin_index) ~= nb_pixel)
                
                axon_collection(k).data = current_myelin;
               
                current_g_ratio = sqrt(length(current_axon_index) / (length(current_myelin_index) + length(current_axon_index)));
                
                axon_collection(k).gRatio = current_g_ratio;
                g_list(k) = current_g_ratio;
                
                mean_g_ratio = mean(g_list);
                if (mean_g_ratio < expected_g_ratio)
                    stop = 1;
                end
            else
                pdf_list(k) = 0;
            end
                
        case 'shrink'
%             map_bound = zeros(dims);
            
            bound = bwboundaries(myelin_map, 'noholes');
            bound = bound{1};
            bound_index = sub2ind(dims,bound(:,1),bound(:,2));
            bound_index = sort(unique(bound_index));
%             map_bound(bound_index) = 1;
            
            axon_map = zeros(dims);
            axon_map(old_axon_index) = 1;
            
            dilate_axon_map = imdilate(axon_map,  ones(3,3));
            smoothed_axon_index = find(dilate_axon_map);
            
            current_myelin_index = setdiff(old_myelin_index, smoothed_axon_index);
            current_myelin_index = union(current_myelin_index, bound_index);
            
            [current_myelin(:,1), current_myelin(:,2)] = ind2sub(dims,current_myelin_index);
            current_axon = myelin2axon(current_myelin);
            current_axon_index = sub2ind(dims,current_axon(:,1),current_axon(:,2));
            
            axon_collection(k).data = current_myelin;
            
            current_g_ratio = sqrt(length(current_axon_index) / (length(current_myelin_index) + length(current_axon_index)));
            
            axon_collection(k).gRatio = current_g_ratio;
            g_list(k) = current_g_ratio;
            
            if (length(current_myelin_index) == length(bound_index))
                pdf_list(k) = 0;
            end
            
            mean_g_ratio = mean(g_list);
            if (mean_g_ratio > expected_g_ratio)
                stop = 1;
            end
    end  
end
end
