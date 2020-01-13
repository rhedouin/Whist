function axonCollection = convertAxonDataToRoundValues(axonCollection)

for k = 1:length(axonCollection)
    axonCollection(k).Centroid = round(axonCollection(k).Centroid);
    axonCollection(k).data = round(axonCollection(k).data);
end
end