function rot = createRotationMatrixAroundAxis(vec, theta)
    % 3D rotation matrix around an axis vec with a theta angle
    % Formula comming from 
    % https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    eps = 1e-4;
    if ((norm(vec) - 1) > eps)
        error('the input should be an unit vector (norm ~= 1)')
    end    
        
    rot = [cos(theta) + (1-cos(theta))*vec(1)^2, vec(1)*vec(2)*(1-cos(theta)) - vec(3)*sin(theta), vec(1)*vec(3)*(1-cos(theta)) + vec(2)*sin(theta); ...
           vec(1)*vec(2)*(1-cos(theta)) + vec(3)*sin(theta), cos(theta) + (1-cos(theta))*vec(2)^2, vec(2)*vec(3)*(1-cos(theta)) - vec(1)*sin(theta); ...
           vec(1)*vec(3)*(1-cos(theta)) - vec(2)*sin(theta), vec(2)*vec(3)*(1-cos(theta)) + vec(1)*sin(theta), cos(theta) + (1-cos(theta))*vec(3)^2];

end
