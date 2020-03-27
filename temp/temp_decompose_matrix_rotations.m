clear
close all

load('theorically_good_10_rotations.mat')

k = 10;

[X, Y, Z, angle_x, angle_y, angle_z]  = decomposeRotationMatrix(rotations(:,:,k));

x_axis_rotated = rotations(:,:,k) * [1; 0; 0];
y_axis_rotated = rotations(:,:,k) * [0; 1; 0];
z_axis_rotated = rotations(:,:,k) * [0; 0; 1];

angle_x = angle_x*(2*15/pi)
angle_y = angle_y*(2*15/pi)
angle_z = angle_z*(2*15/pi)

visualize3dArrow(x_axis_rotated', 0, 'b');
visualize3dArrow(y_axis_rotated', 0, 'k');
visualize3dArrow(z_axis_rotated', 0, 'r');

% end