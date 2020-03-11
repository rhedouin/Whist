clear
close all

base_folder = ['/project/3015069.04/'];
rotations_folder = [base_folder 'data/rotations/'];
 
load([rotations_folder '20_fiber_orientations_no_rotation.mat']);

rotation_axis = ParticleSampleSphere('N',16);
visualize3dArrow(rotation_axis)

unit = [0; 0; 1];
nb_rotations = 16;
for l = 1:20
    for k = 1:nb_rotations
        rotations(:,:,k) = createRotationMatrixFromVector(rotation_axis(k,:));

        alldirections(k,l,:) = rotations(:,:,k) * fiber_directions(l,:)';
        theta(l,k) = acos(sign(alldirections(k,l,3))*alldirections(k,l,3));

    end    
end

figure
for k = 1:nb_rotations
    theta_sort(:,k) = sort(theta(:,k));
    plot(theta_sort(:,k) )
    hold on
end

save([rotations_folder 'theorically_good_' num2str(nb_rotations) '_rotations.mat'], 'rotations')


return



















unit = [0;0;1];
% 
diff = -pi/2;
it = 0

for j = 1:9
    [X(:,:,j) Y(:,:,j) Z(:,:,j) angle_x(j) angle_y(j) angle_z(j)] = decomposeRotationMatrix(best_rotations(:,:,j));
    
    orientations(j,:) = best_rotations(:,:,j) * unit;
%     orientations_without_z(j,:) = Y(:,:,j)*X(:,:,j)*unit;
end


angle_x_deg = rad2deg(angle_x)
angle_y_deg = rad2deg(angle_y)
angle_z_deg = rad2deg(angle_z)

angle_x_cm = angle_x*16/(pi/2)
angle_y_cm = angle_y*16/(pi/2)
angle_z_cm = angle_z*16/(pi/2)

visualize3dArrow(orientations);

for l = 1:20
    for j = 1:9
        current_dir = best_rotations(:,:,j) * rotations20(:,:,l) * unit;
        theta(l,j) = acos(sign(current_dir(3))*current_dir(3));
        
        current_dir_without_z = Y(:,:,j) * X(:,:,j) * rotations20(:,:,l) * unit;
        theta_without_z(l,j) = acos(sign(current_dir_without_z(3))*current_dir_without_z(3));
             
    end
    theta_sort(l,:) = sort(theta(l,:));
    figure(2)
    plot(theta_sort(l,:))
    hold on;
    title('theta sort')

end



for j = 1:9
    [X(:,:,j) Y(:,:,j) Z(:,:,j) angle_x(j) angle_y(j) angle_z(j)] = decomposeRotationMatrix(new_rotations(:,:,j));
    
    orientations(j,:) = new_rotations(:,:,j) * unit;
end


angle_x_deg = rad2deg(angle_x)
angle_y_deg = rad2deg(angle_y)
angle_z_deg = rad2deg(angle_z)

angle_x_cm = angle_x*16/(pi/2)
angle_y_cm = angle_y*16/(pi/2)
angle_z_cm = angle_z*16/(pi/2)

visualize3dArrow(orientations);



for l = 1:20
    for j = 1:9
        current_dir = new_rotations(:,:,j) * rotations20(:,:,l) * unit;
        theta(l,j) = acos(sign(current_dir(3))*current_dir(3));
        
        current_dir_without_z = Y(:,:,j) * X(:,:,j) * rotations20(:,:,l) * unit;
        theta_without_z(l,j) = acos(sign(current_dir_without_z(3))*current_dir_without_z(3));
             
    end
    theta_sort(l,:) = sort(theta(l,:));
    figure(4)
    plot(theta_sort(l,:))
    hold on;
    title('theta sort')

end

keyboard;














for j = 1:9
    test(j,:) = best_rotations(:,:,j) * unit;
    test_without_z(j,:) = Y(:,:,j) * X(:,:,j) * unit;
end
 
visualize3dArrow(test);
visualize3dArrow(test_without_z);
visualize3dArrow(best_orientations);

for j=1:9
    toto(j) = acos(abs(test_without_z(j,2)))
end
for j=1:9
    toto2(j) = acos(abs(test(j,2)))
end
return







keyboard;

for j = 1:9
    exvivo_orientations(j,:) = rotations9(:,:,j) *unit;
    exvivo_axis_rotations(j,:) = vrrotmat2vec(rotations9(:,:,j));
                 
    best_axis_rotations(j,:) = vrrotmat2vec(best_rotations(:,:,j));
    classic_rotations(:,:,j) = createRotationMatrixFromVector(best_orientations(j,:));
    classic_axis_rotations(j,:) = vrrotmat2vec(classic_rotations(:,:,j));
end

visualize3dArrow(exvivo_axis_rotations(:,1:3));
visualize3dArrow(best_orientations);
visualize3dArrow(best_axis_rotations(:,1:3));
visualize3dArrow(classic_axis_rotations(:,1:3));
visualize3dArrow(classic_axis_rotations(:,4) .* classic_axis_rotations(:,1:3));

return

% 
% figure
% for l = 1:20
%     for j = 1:9
%         current_dir = new_rotations(:,:,j) * rotations20(:,:,l) * unit;
%         theta(l,j) = acos(sign(current_dir(3))*current_dir(3));
%         
%     end
%     theta_sort(l,:) = sort(theta(l,:));
%     plot(theta_sort(l,:))
%     hold on
%     title('9 best rotations')
% end
% 
% 
% figure
% for l = 1:20
%     for j = 1:9
%         current_dir = new_rotations(:,:,j) * rotations20(:,:,l) * unit;
%         theta(l,j) = acos(sign(current_dir(3))*current_dir(3));
%         
%     end
%     theta_sort(l,:) = sort(theta(l,:));
%     plot(theta_sort(l,:))
%     hold on
% end
% % visualize3dArrow(current_dir);
% 
% 
% max_first = max(theta_sort(:,1));
% min_last = min(theta_sort(:,end));
% new_diff = min_last - max_first
% 
% figure
% for l = 1:20
%     for j = 1:9
%         current_dir = rotations9(:,:,j) * rotations20(:,:,l) * unit;
%         directions3(j,:) = current_dir;
%         theta(l,j) = acos(sign(current_dir(3))*current_dir(3));
%     end
%     theta_sort(l,:) = sort(theta(l,:));
%     
%     max_first = max(theta_sort(:,1));
%     min_last = min(theta_sort(:,end));
%     new_diff = min_last - max_first;
%     
%     plot(theta_sort(l,:))
%     hold on
%     title('ex-vivo')
% end
% max_first = max(theta_sort(:,1));
% min_last = min(theta_sort(:,end));
% new_diff = min_last - max_first
% return

for nb = 1:500

    nb_top_orientations = 0;
    while nb_top_orientations ~= 9
        random_orientations = ParticleSampleSphere('N',18,'Nitr',100);
        positive_values = random_orientations(:,3);
        random_orientations = random_orientations(find(positive_values > 0),:);
        nb_top_orientations = size(random_orientations,1);
    end
    
    for j = 1:9
        random_rotations(:,:,j) = RotMatrix(pi/2, random_orientations(j,:));
    end
    
    for l = 1:20
        for j = 1:9
            current_dir = random_rotations(:,:,j) * rotations20(:,:,l) * unit;
            theta(l,j) = acos(sign(current_dir(3))*current_dir(3));
        end
        theta_sort(l,:) = sort(theta(l,:));
    end
    
    max_first = max(theta_sort(:,1));
    min_last = min(theta_sort(:,end));
    new_diff = min_last - max_first
    if new_diff > diff
        it = it+1
        diff = new_diff
        best_rotations = random_rotations;
        best_orientations = random_orientations;
        best_theta_sort = theta_sort;
    end
end
figure
for j = 1:20
    plot(best_theta_sort(j,:))
    hold on
end
ylabel('theta')
xlabel('rotations')
title('rotation optimized, diff = 0.86')
visualize3dArrow(best_orientations)

save('/project/3015069.01/model/code/mathTools/9rotations_optimized.mat', 'best_rotations', 'best_orientations')
return
    
% for j = 1:9
%     current_dir = rotations9(:,:,j) * unit;
%     directions(j,:) = current_dir;
% end

% for j = 1:9
%     current_dir = rotations20(:,:,15) * rotations9(:,:,j) * unit;
%     directions2(j,:) = current_dir;
% end





visualize3dArrow(directions3)
title('ex-vivo')

figure
for l = 1:20
    for j = 1:9
        current_dir = rotations9(:,:,j) * rotations20(:,:,l) *  unit;
        directions3(j,:) = current_dir;
        theta(l,j) = acos(sign(current_dir(3))*current_dir(3));
    end
    theta_sort(l,:) = sort(theta(l,:))
    plot(theta_sort(l,:))
    hold on
    title('current rotation, diff 0.08')

end
visualize3dArrow(directions3)











orientation_all_sphere = ParticleSampleSphere('N',10);
for j = 1:10  
    rotation10_all_sphere(:,:,j) = createRotationMatrixFromVector(orientation_all_sphere(j,:));
end
visualize3dArrow(orientation_all_sphere)

figure

clear theta_sort theta directions3
for l = 1:1
    for j = 1:10
        current_dir = rotation10_all_sphere(:,:,j) * rotations20(:,:,l) *  unit;
        directions3(j,:) = current_dir;
        theta(l,j) = acos(sign(current_dir(3))*current_dir(3));
    end
    theta_sort(l,:) = sort(theta(l,:))
    plot(theta_sort(l,:))
    hold on
    title('all_sphere')

end
visualize3dArrow(directions3)
title('all_sphere')


for j = 1:10
    test(j,:) = rotation10_all_sphere(:,:,j)*unit;
end
visualize3dArrow(test)

keyboard;

for j = 1:9
    current_dir = rotations9(:,:,j) *  unit;
    directions3(j,:) = current_dir;
end


sum(directions.*directions,2)
sum(directions.*directions2,2)
sum(directions.*directions3,2)

for j = 1:9
    rotations_new(:,:,j) = createRotationMatrixFromVector(rotations(j,:))
end

