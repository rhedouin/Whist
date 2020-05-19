clear
close all

load('theorically_good_10_rotations.mat')

for k = 1:10
    
     [~, ~, ~, angle_x(k), angle_y(k), angle_z(k)] = decomposeRotationMatrix(rotations(:,:,k));
     x_rotated = rotations(:,:,k)* [1; 0; 0];
     y_rotated = rotations(:,:,k)* [0; 1; 0];
     z_rotated = rotations(:,:,k)* [0; 0; 1];
     
     new_spot = [x_rotated'; y_rotated'; z_rotated'];
     
     visualize3dArrow(new_spot, ')

end