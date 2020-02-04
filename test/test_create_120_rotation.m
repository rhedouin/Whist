% Create rotation matrix 
clear
close all
vec = randn(3,1);
vec = vec / norm(vec)

ang = 2*pi/3;
rot = rotationmat3D(ang, vec)

current_vec = [0; 0; 1];

for  k = 1:3
    save_dir(k, :) = current_vec;
    current_vec = rot * current_vec;
end

visualize3dArrow(save_dir)
% vec should be the baricenter
visualize3dArrow(vec', 0, 'r')