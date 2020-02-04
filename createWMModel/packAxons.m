function [axon_collection, FVF] = packAxons(axon_collection, mask, max_iteration, max_FVF, step)
% author : Renaud Hedouin 
% email : renaud.hedouin@gmail.com
% Code adapted from Tom Mingasson toolbox 

dims = size(mask);

disp(' ')
disp('Packing in process...')
disp(' ')

N = length(axon_collection);
FVF(1) = 0;

plot = 1;

for iter = 1:max_iteration
    disp(['Iteration : ' num2str(iter)]);
    if (mod(iter, 10)) == 0
        [~, ~, FVF, ~] = createModelFromData(axon_collection, mask, plot);
        
        disp(['FVF : ' num2str(FVF)]);
    
        pause(0.1)
        if(FVF > max_FVF)
            break;
        end
    end
    axon_collection = packAxonsOneIteration(axon_collection, dims, FVF, step);

end
end


function axon_collection = packAxonsOneIteration(axon_collection, dims, FVF, step)

pts = cat(1,axon_collection(:).Centroid);
N=size(pts,1);

% Progressive attraction coefficient
Kcenter0 = step*(1 - FVF);     % center step coeff for disks withOUT overlapping
Kcenter1 = 0;         % center step coeff for disks with overlapping
Krep = step*(1-FVF)/2;

% check all pair of axon intersections
Lbin = computeAxonsSuperpositionMatrix(axon_collection, dims);
Lbin2 = repmat(Lbin, [1, 1, 2]);
inter1_index = repmat(sum(Lbin,2),[1,2]);
inter1_index(inter1_index>0)=1;   % disks that overlap
inter0_index = 1 - inter1_index;  % disks that NOT overlap

% attraction
pts_centered = dims/2-pts;
attraction_norm = sqrt(pts_centered(:,1).^2+pts_centered(:,2).^2);
attraction = pts_centered./repmat(attraction_norm,[1,2]);

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

MyGrad = Kcenter0.*attraction.*inter0_index + Kcenter1.*attraction.*inter1_index + Krep.*repulsion.*inter1_index ;

for k = 1:N
    axon_collection(k).Centroid = axon_collection(k).Centroid + MyGrad(k,:);
    axon_collection(k).data = axon_collection(k).data + MyGrad(k,:);
end
end



