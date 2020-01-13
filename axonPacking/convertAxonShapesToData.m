function axonCollection = convertAxonShapesToData(axon_shapes, pts)

for k = 1:length(axon_shapes)
    axonCollection(k).data = axon_shapes{k} + pts(k,:);
end
end