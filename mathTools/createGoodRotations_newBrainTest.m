clear
close all

max_theta_diff = 0;

base_folder = ['/project/3015069.04/'];
rotations_folder = [base_folder 'data/rotations/'];

load([rotations_folder '20_fiber_orientations_rotationnot_spread_evenly.mat']);
load([rotations_folder 'theorically_good_10_rotations.mat']);

nb_rotations = 10;

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
theta_diff = min(theta_sort(:, end)) - max(theta_sort(:, 1))
























