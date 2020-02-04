% Ex vivo post processing 
% Concatenate signals
base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/';
dti_folder = [base_folder 'dwi/orientation-9/dti/'];
input_gre_folder = [base_folder 'gre_kwok/'];
output_gre_folder = [base_folder 'gre_renaud/'];

fa_list = {'05','10','15','20','35','60'};
lfa = length(fa_list);

for korient = 1:9
    orient = ['orientation-' num2str(korient)]
    input_orientation_folder = [input_gre_folder orient '/'];
    output_orientation_folder = [output_gre_folder orient '/'];
    mkdir(output_orientation_folder);
  
    input_v1_path = [dti_folder 'BrainSample-2_ses-03_dwi_orientation-9_NSA-16_log-euclidean_dti_V1_2_BrainSample-2_ses-03_gre_' orient '_fa-20_magn_echo-01.nii.gz'];
    
    v1_nii = load_untouch_nii(input_v1_path);
    
    v1_img = v1_nii.img;
    v1_img(v1_img > 1) = 1;
    v1_img(v1_img < -1) = -1;
    
    theta = acos(abs(v1_img(:,:,:,3)));
    
    input_mask_path = [input_orientation_folder 'BrainSample-2_ses-03_gre_' orient '_brain_mask.nii.gz'];
    mask_nii = load_untouch_nii(input_mask_path);
    mask = single(mask_nii.img);
    dims = size(mask);
    
    theta = mask .* theta;
    
    mask_nii.img = theta;
    mask_nii.hdr.dime.datatype = 16;
    mask_nii.hdr.dime.bitpix = 32;
      
    output_theta_path = [output_orientation_folder 'BrainSample-2_ses-03_gre_' orient '_theta.nii.gz'];
    save_untouch_nii(mask_nii, output_theta_path);

    for kfa = 1:lfa
        tic()
        fa = ['fa-' fa_list{kfa}]
        input_fa_folder = [input_orientation_folder fa '/'];
        output_fa_folder = [output_orientation_folder fa '/'];
      
        mkdir(output_fa_folder);
        
        input_magn_path = [input_fa_folder 'BrainSample-2_ses-03_gre_' orient '_' fa '_magn.nii.gz'];
        input_phase_path = [input_fa_folder 'local-field-lbv/BrainSample-2_ses-03_gre_' orient '_' fa '_phase.nii.gz'];
        
        magn_nii = load_untouch_nii(input_magn_path);
        magn = single(magn_nii.img);
        phase = single(load_nii_img_only(input_phase_path));
        
        load([input_fa_folder 'BrainSample-2_ses-03_gre_' orient '_' fa '_header.mat'])
        
        nb_TE = length(TE);
        
        %% phase unwrapping
        [~,~,unwarped_phase] = estimateTotalField(phase, magn, matrixSize,voxelSize,...
            'Method','Optimum weights','Unwrap','bestpath3d',...
            'unit','Hz','TE',TE,'B0',B0,...
            'mask',mask);
        
        unwarped_phase = unwarped_phase.*repmat(mask, [1 1 1 nb_TE]);
        
        %% magnitude and phase normalization
        normed_magn = zeros([dims nb_TE]);
        normed_phase = zeros([dims nb_TE-2]);

        for k = 1:dims(1)
            for l = 1:dims(2)
                for m = 1:dims(3)
                    if (mask(k,l,m) ~= 0)
                        
                        normed_magn(k,l,m,:) = abs(magn(k,l,m,:)) / abs(magn(k,l,m,1));
                        
                        temp_phase = squeeze(unwarped_phase(k,l,m,:))';
                        normed_phase(k,l,m,:) = temp_phase(3:end) - temp_phase(1) - (temp_phase(2)-temp_phase(1))*(TE(3:end) - TE(1))/(TE(2) - TE(1));
                    end
                end
            end
        end
        

        concatenate_signal = cat(4, theta, normed_magn, normed_phase);
        
        %% Save all images
        output_unwarped_phase_path = [output_fa_folder 'BrainSample-2_ses-03_' orient '_' fa '_unwrapped-phase-masked.nii.gz'];
        output_normed_phase_path = [output_fa_folder 'BrainSample-2_ses-03_' orient '_' fa '_normed-phase-masked.nii.gz'];
        output_normed_magn_path = [output_fa_folder 'BrainSample-2_ses-03_' orient '_' fa '_normed-magn-masked.nii.gz'];
        output_concatenate_signal_path = [output_fa_folder 'BrainSample-2_ses-03_' orient '_' fa '_concatenate_signal.nii.gz'];
        
        magn_nii.hdr.dime.datatype = 16;
        magn_nii.hdr.dime.bitpix = 32;
        magn_nii.hdr.dime.dim(1) = 4;
        magn_nii.hdr.dime.dim(5) = 12;
        magn_nii.img = single(unwarped_phase);
        save_untouch_nii(magn_nii, output_unwarped_phase_path);
     
        magn_nii.hdr.dime.dim(5) = 12;
        magn_nii.img = single(normed_magn);
        save_untouch_nii(magn_nii, output_normed_magn_path);
        
        magn_nii.hdr.dime.dim(5) = 10;
        magn_nii.img = single(normed_phase);
        save_untouch_nii(magn_nii, output_normed_phase_path);

        magn_nii.hdr.dime.dim(5) = 23;
        magn_nii.img = single(concatenate_signal);
        save_untouch_nii(magn_nii, output_concatenate_signal_path);
        toc()
        
    end   
end

    