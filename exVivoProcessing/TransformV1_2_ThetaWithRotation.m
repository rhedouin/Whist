% Transform Magn Phase in Polyfit cartesian
function theta =  TransformV1_2_ThetaWithRotation(V1, rotation, mask)

dims = size(V1);
theta = zeros(dims(1:3));

eps = 1e-4;

for k = 1:dims(1)
    for l = 1:dims(2)
        for m = 1:dims(3)
            if mask(k, l, m)
                current_V1 = squeeze(V1(k,l,m,:));
                if sum(current_V1) > eps
                    V1_rotated = rotation * current_V1;
                else
                    V1_rotated = rotation * [0;0;1];
                end
                
                if V1_rotated(3) > 1
                    V1_rotated(3) = 1;
                elseif V1_rotated(3) < -1
                    V1_rotated(3) = -1;
                end
                
                theta(k,l,m) = acos(abs(V1_rotated(3)));
                              
            end
        end
    end
end
end