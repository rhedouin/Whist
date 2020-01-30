function [axon_shapes, pts] = convertAxonDataToShapes(axonCollection)

for k = 1:length(axonCollection)
    pts(k,:) = axons(k).Centroid;
    axon_shapes{k} = axons(k).data - pts(k,:);
end
end