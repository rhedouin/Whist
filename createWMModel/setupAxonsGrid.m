function  [axonCollection, dims] = setupAxonsGrid(axonCollection, dims)

maxSize = 0;

N = length(axonCollection);
for k = 1:N
    currentSize = max(max(axonCollection(k).data)  - min(axonCollection(k).data ));
    if (currentSize > maxSize)
        maxSize = max(currentSize);
    end
end
maxRadius = (maxSize /2) + 5;
% dimension of the square area
sqrt_N = round(sqrt(N));

% Random positions on a grid for the N axons
if ~exist('dims')
    space = round((16/10)*maxRadius); 
    side = 2*maxRadius + (sqrt_N-1) * space;
    grid_spacing = maxRadius : space : side;

    dims = [side side];
else
    grid_spacing = round(linspace(maxRadius, dims(1) - maxRadius, sqrt_N));
end
    
[X_grid Y_grid] =  meshgrid(grid_spacing, grid_spacing);

X_grid =X_grid(:);
Y_grid =Y_grid(:);
Permutations = randperm(N);
pts = zeros(N,2);

for k=1:N
    pts(k,:) = [X_grid(Permutations(k)) Y_grid(Permutations(k))];
    axonCollection(k).data = axonCollection(k).data - axonCollection(k).Centroid + pts(k,:);
    axonCollection(k).Centroid =  pts(k,:);
end
end




