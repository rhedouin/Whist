function [orientations, dispersion] = generateVMFsample_for_dispersion(avrg_orient, kappa, expected_dispersion, nb_orientations)
%generate vMF sample 

rot = createRotationMatrixFromVector(avrg_orient);

if kappa > 1000
    kappa = 1000;
end

tol = 0.01;
iteration = 0;
while true 
    
    iteration = iteration + 1;
    if (iteration > 10000)
        error('Too much iteration');
    end
    sum_dot_prod = 0;
    
    clear orientations
    
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
        
        sum_dot_prod = sum_dot_prod + dot(orientations(l,:), avrg_orient).^2;
        
    end
    
    rho(iteration) = 1 - (sum_dot_prod / nb_orientations);
    
    if (abs(rho(iteration) - expected_dispersion) < tol)
        dispersion = rho(iteration)
        break;
    end
      
end

end



