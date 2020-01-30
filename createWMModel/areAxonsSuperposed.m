function out = areAxonsSuperposed(axon1, axon2, dims)

a1min = min(axon1);
a1max = max(axon1);

a2min = min(axon2);
a2max = max(axon2);

if ((sum(a1min > a2max) >= 1) || (sum(a1max < a2min) >= 1))
    out = 0;
    return;
end

ind1 = round(axon1);
sub1 = sub2ind(dims,ind1(:,1),ind1(:,2));  

ind2 = round(axon2);
sub2 = sub2ind(dims,ind2(:,1),ind2(:,2));  

C = intersect(sub1,sub2);
if (length(C) > 0)
    out = 1;
else 
    out = 0;
end
end


