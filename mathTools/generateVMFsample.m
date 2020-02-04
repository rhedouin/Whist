function orientations = generateVMFsample(avrg_orient, kappa, nb_orientations)
%generate vMF sample 

rot = createRotationMatrixFromVector(avrg_orient);

if kappa > 1000
    kappa = 1000;
end

for l = 1:nb_orientations
    y = rand(1);
    
    if kappa < 50
        ck = (exp(kappa) - exp(-kappa))/kappa;
        w = (1/kappa)*log(exp(-kappa) + kappa*y*ck);
    else
        w = (log(y) + kappa)/kappa;
    end
    
    theta = 2*pi*rand(1);
    vec = [cos(theta) sin(theta)];
    
    v(l,:) = [sqrt(1-w^2)*vec, w];
    orientations(l,:) = rot * v(l,:)';
    
end
end



