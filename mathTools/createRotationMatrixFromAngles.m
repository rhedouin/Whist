function rot = createRotationMatrixFromAngles(theta, phi)
    % Already checked (should works)

%     rot_y = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
    rot_y = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    rot_z = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

    rot = rot_z * rot_y;

end