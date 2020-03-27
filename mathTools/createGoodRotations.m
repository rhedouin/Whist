clear
close all

max_theta_diff = 0;

base_folder = ['/project/3015069.04/'];
rotations_folder = [base_folder 'data/rotations/'];

load([rotations_folder '20_fiber_orientations_no_rotation.mat']);

nb_rotations = 16;

for num = 1:1000
    
    rotation_axis = ParticleSampleSphere('N',nb_rotations);
    % visualize3dArrow(rotation_axis)
    
    unit = [0; 0; 1];
    nb_fiber_orientations = 20;
    
    for l = 1:nb_fiber_orientations
        for k = 1:nb_rotations
%             rotations(:,:,k) = createRotationMatrixAroundAxis(rotation_axis(k,:), pi/2);
            
            alldirections(k,l,:) = rotations(:,:,k) * fiber_directions(l,:)';
            theta(l,k) = acos(sign(alldirections(k,l,3))*alldirections(k,l,3));     
        end
    end
    
    figure
    
    for k = 1:nb_fiber_orientations
        theta_sort(k, :) = sort(theta(k, :));
        plot(theta_sort(k, :) )
        hold on
    end
    num
    theta_diff = min(theta_sort(:, end)) - max(theta_sort(:, 1))
    
    if theta_diff > max_theta_diff
        'toto'
        best_rotations = rotations;
        max_theta_diff = theta_diff;
    end
end
rotations = best_rotations;
save([rotations_folder 'theorically_good_' num2str(nb_rotations) '.mat'], 'rotations')

return

























