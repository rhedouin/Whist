function axon_collection = avoidAxonOverlap(axon_collection, dims)

Lbin = computeSuperpositionMatrix(axon_collection, dims);
nb_axons_overlap = sum(Lbin, 'all')/2;

iter = 0;

while and((nb_axons_overlap ~= 0), iter < 50);
    
    iter = iter + 1;
    display(['Number of axons which overlap : ' num2str(nb_axons_overlap)]);
    pts = cat(1,axon_collection(:).Centroid);
    % Test progressive attraction coefficient
        
    Krep = 0.5;
    
    N=size(pts,1);
    % intersection
    Lbin = computeSuperpositionMatrix(axon_collection, dims);
    nb_axons_overlap = sum(Lbin, 'all')/2;

    Lbin2 = repmat(Lbin, [1, 1, 2]);
    inter1_index = repmat(sum(Lbin,2),[1,2]);
    inter1_index(inter1_index>0)=1;   % disks that overlap
    inter0_index = 1 - inter1_index;  % disks that NOT overlap
    
    % repulsion
    
    pts_replic = permute(repmat(pts,[1,1,N]), [1 3 2]);
    pts_replic_switch = permute(pts_replic, [2 1 3]);
    
    U = pts_replic_switch - pts_replic;
    
    Usum  = squeeze(sum(U.*Lbin2, 1));
    Unorm = sqrt(Usum(:,1).^2 + Usum(:,2).^2);
    Unorm(Unorm==0) = 1;
    Unormalization = repmat(Unorm,[1,2]);
    
    Usum_normed = Usum./Unormalization;
    repulsion = Usum_normed;
    
    MyGrad = Krep.*repulsion.*inter1_index ;
    
    for k = 1:N
        axon_collection(k).Centroid = axon_collection(k).Centroid + MyGrad(k,:);
        axon_collection(k).data = axon_collection(k).data + MyGrad(k,:);
    end
end
end










