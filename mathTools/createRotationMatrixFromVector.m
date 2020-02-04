function [rot, theta, phi] = createRotationMatrixFromVector(vec)
    % Already checked (should works)
    % This function take a 3d vector as an input and give as output the
    % rotation matrix that send [0;0;1] on this vector 
    % Spherical coordinate with physic convention (see wiki)  
    eps = 1e-4;
    if ((norm(vec) - 1) > eps)
        error('the input should be an unit vector (norm ~= 1)')
    end
    
    if (dot(vec, [0 0 1]) == 1)
        rot = eye(3);
        theta = 0;
        phi = 0;
        return;
    end

    unit = vec / norm(vec);
    
    phi = atan(unit(2)/unit(1));
    theta = sign(unit(1))*acos(unit(3));
        
    rot_y = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    rot_z = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

    rot = rot_z*rot_y;

end
