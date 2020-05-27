% Ex vivo post processing 
% Concatenate signals
input_folder = '/project/3015069.01/derived/BS-3/ses-03/';
dti_folder = [input_folder 'dwi/orientation-9/dti/'];

fa_list = [5, 10, 20, 30, 40, 50, 70];
lfa = length(fa_list);

for kOrient = 1:10
    gre_orientation_folder = [input_folder 'gre_renaud/' orient '/'];
    
    input_v1_path = [dti_folder 'BrainSample-2_ses-03_dwi_orientation-9_NSA-16_log-euclidean_dti_V1_2_BrainSample-2_ses-03_gre_' orient '_fa-20_magn_echo-01.nii.gz'];
    
    v1_nii = load_untouch_nii(input_v1_path);
    
    v1_img = v1_nii.img;
    v1_img(v1_img > 1) = 1;
    v1_img(v1_img < -1) = -1;
    
    theta = acos(abs(v1_img(:,:,:,3)));
    
    for kfa = 1:lfa
        
        fa = ['fa-' fa_list{kfa}]
        fa_folder = [gre_orientation_folder fa '/'];
        
        input_magn_path = [fa_folder 'BrainSample-2_ses-03_gre_' orient '_' fa '_normed-magn-masked.nii.gz'];
        input_norm_path = [fa_folder 'BrainSample-2_ses-03_gre_' orient '_' fa '_normed-phase-masked.nii.gz'];
        
        magn = load_nii_img_only(input_magn_path);
        norm = load_nii_img_only(input_norm_path);
        concatenate_signal = cat(4, theta, magn, norm);
        
        v1_nii.hdr.dime.dim(1) = 4;
        v1_nii.hdr.dime.dim(5) = 23;
        v1_nii.img = concatenate_signal;
        
        output_path = [fa_folder 'BrainSample-2_ses-03_' orient '_' fa '_concatenate_signal.nii.gz'];
        save_untouch_nii(v1_nii, output_path);
        
    end   
end

    