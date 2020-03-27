clear
close all

load('theorically_good_10_rotations.mat')

angle_x = pi/6;
angle_y = pi/3;
angle_z = -pi/4;

X = [1, 0, 0; 0, cos(angle_x), -sin(angle_x); 0, sin(angle_x), cos(angle_x)];
Y = [cos(angle_y), 0, sin(angle_y); 0, 1, 0; -sin(angle_y), 0 cos(angle_y)];
Z = [cos(angle_z), -sin(angle_z), 0; sin(angle_z), cos(angle_z), 0; 0, 0, 1];

rotations = Z*Y*X;

[X, Y, Z, angle_x, angle_y, angle_z]  = decomposeRotationMatrix(rotations);

x_axis_rotated = rotations * [1; 0; 0];
y_axis_rotated = rotations * [0; 1; 0];
z_axis_rotated = rotations * [0; 0; 1];

angle_x = angle_x*(30/pi)
angle_y = angle_y*(30/pi)
angle_z = angle_z*(30/pi)

visualize3dArrow(x_axis_rotated', 0, 'b')
visualize3dArrow(y_axis_rotated', 0, 'k')
visualize3dArrow(z_axis_rotated', 0, 'r')

% end