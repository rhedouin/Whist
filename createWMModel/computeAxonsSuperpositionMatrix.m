function Dist = computeAxonsSuperpositionMatrix(axon_collection, dims)
matrix_dims = length(axon_collection);
Dist = zeros(matrix_dims, matrix_dims);
for k = 1:matrix_dims
    for l = k+1:matrix_dims
        Dist(k,l) = areAxonsSuperposed(axon_collection(k).data, axon_collection(l).data, dims);
        Dist(l,k) = Dist(k,l);
    end
end
end