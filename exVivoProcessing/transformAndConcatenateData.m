% Transform and concatenate BS3 gre dwi data 
function out = transformAndConcatenateBS3Data(fa)

base_folder = '/project/3015069.01/derived/BS-3/';
gre_folder = [base_folder 'ses-mri03/gre/lowres/'];
gre_ref_folder = [gre_folder 'ref_space/'];
concatenate_folder = [gre_ref_folder 'concatenate_signals/'];

% Here you need to give the rotations corresponding to your different acquisition orientations
load('/project/3015069.04/data/rotations/BS3_highres_rotations_ref_2_orientations.mat')

% Here goes your dti main orientations (V1 image) to extract the theta angle with the B0 orientation
V1 = load_nii_img_only('/project/3015069.01/derived/BS-3/ses-mri02/dwi/results/BS-3_ses-mri02_dwi_lpca_eddy_NSA-20_dti_V1_2_gre_lowres_ref.nii.gz');

nb_TE = 18;
nb_orientations = 10;

% Here goes your mask
mask = load_nii_img_only([concatenate_folder 'BS-3_ses-mri03_acq-lowres_gre_gradunwarp_mask_sphere_ref.nii.gz']);

    concatenate_signals = [];
    for kOri = 1:nb_orientations
        orient_folder = [gre_ref_folder 'orientation-' num2str(kOri) '/'];
        mask_orient_folder = [gre_folder 'orientation-' num2str(kOri) '/'];
            
        fa_folder = [orient_folder 'fa-' num2str(fa) '/'];
        
% Here goes your magnitude and phase images 
        magn_nii = load_untouch_nii([fa_folder 'BS-3_ses-mri03_acq-lowres_orientation-' num2str(kOri) '_fa-' num2str(fa) '_part-mag_gre_gradunwarp_2_ref.nii.gz']);
        unwrapped_phase_nii = load_untouch_nii([fa_folder 'BS-3_ses-mri03_acq-lowres_orientation-' num2str(kOri) '_fa-' num2str(fa) '_gre_gradunwarp_unwrapped-phase_2_ref.nii.gz']);

        magn = single(magn_nii.img);
        unwrapped_phase = single(unwrapped_phase_nii.img);
        
        display('transform V1 to theta with rotation')
        tic()
        theta =  TransformV1_2_ThetaWithRotation(V1, rotations(:,:,kOri), mask);
        toc()
     
        display('transform magn and unwrapped phase to cartesian polyfit')
        tic()
        [real_complex_polyfit, imag_complex_polyfit] = TransformMagnPhase_2_PolyfitCartesian(magn, unwrapped_phase, mask);
        toc()
        
        concatenate_signals = cat(4,concatenate_signals, theta, real_complex_polyfit, imag_complex_polyfit);
                
    end
    magn_nii.hdr.dime.dim(5) = (2*nb_TE + 1) * nb_orientations;
    magn_nii.img = single(concatenate_signals);
    magn_nii.hdr.dime.datatype = 16;
    magn_nii.hdr.dime.bitpix = 32;
    
% Here you save your images    
    save_untouch_nii(magn_nii, [concatenate_folder 'BS-3_ses-mri03_acq-lowres_all_orientations_fa-' num2str(fa) '_concatenate_signals_2_ref.nii.gz']);
    
end

        
        
        
