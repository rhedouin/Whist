function [X, Y, Z, angle_x, angle_y, angle_z]   = decomposeRotationMatrix(M)

    % Decompose a general rotation matrix into 3 rotation matrices around
    % the x, y and z axis.
    
    angle_x = atan2(M(3,2), M(3,3));
    angle_y = atan2(-M(3,1), sqrt(M(3,2)*M(3,2) + M(3,3)*M(3,3)));
    angle_z = atan2(M(2,1),M(1,1));
    
    X = [1, 0, 0; 0, cos(angle_x), -sin(angle_x); 0, sin(angle_x), cos(angle_x)];
    Y = [cos(angle_y), 0, sin(angle_y); 0, 1, 0; -sin(angle_y), 0 cos(angle_y)];
    Z = [cos(angle_z), -sin(angle_z), 0; sin(angle_z), cos(angle_z), 0; 0, 0, 1];
    
%     keyboard;
    
%     R = Z*Y*X;    
%     
%     X*[1; 0; 0]
%     X*[0; 1; 0]
%     X*[0; 0; 1]
%     keyboard;

   
end