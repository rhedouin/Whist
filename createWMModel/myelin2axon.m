function axon_index = myelin2axon(myelin_index)

indminx = min(myelin_index(:,1))-1;
indminy = min(myelin_index(:,2))-1;

sizeax = [max(myelin_index(:,1)) - indminx, max(myelin_index(:,2)) - indminy];
% Initialize myelin with small size
myelin_map = zeros(sizeax);
ind = sub2ind(sizeax, myelin_index(:,1) - indminx, myelin_index(:,2) - indminy);
myelin_map(ind) = 1;

% Get axon region (BG region inside myelin mask)
axon_map = imfill(myelin_map, 'holes') - myelin_map;
[xa,ya] = find(axon_map);

if isempty([xa,ya])
    xa=1;
    ya=1;
end

axon_index = zeros(length(xa + indminx),2);
axon_index(:,1) = xa + indminx;
axon_index(:,2) = ya + indminy;

end