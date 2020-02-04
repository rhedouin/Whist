% Ex vivo post processing 
% Concatenate signals


for krot = 2:9

    clear  rotated_V1  theta
    rot = sphere_rotations(:,:,krot);
    dims = [128 128 128]
    for k = 1:128
        k
        for l = 1:128
            for m = 1:128
                rotated_V1(k, l, m, :) = rot*squeeze(original_V1(k,l,m,:));
                
                if rotated_V1(k,l,m,3) > 1
                   rotated_V1(k,l,m,3) = 1
                elseif rotated_V1(k,l,m,3) < -1
                   rotated_V1(k,l,m,3) = -1
                end
                
                theta(k, l, m) = acos(abs(rotated_V1(k, l, m, 3)));
            end
        end
    end
    
    A.hdr.dime.dim(1) = 4;
    A.hdr.dime.dim(5) = 3;
    A.img = rotated_V1;
       
    output_path = ['orientation-' num2str(krot) '/BrainSample-2_ses-03_dwi_orientation-' num2str(krot) '_V1_registered.nii.gz'];
    
    save_untouch_nii(A, output_path);
    
    A.hdr.dime.dim(1) = 3;
    A.hdr.dime.dim(5) = 1;
    A.img = theta;
    
    output_path = ['orientation-' num2str(krot) '/BrainSample-2_ses-03_dwi_orientation-' num2str(krot) '_theta_registered.nii.gz'];

    save_untouch_nii(A, output_path);
end
 
    