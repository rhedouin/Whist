% unwarp phase BS3

base_folder = '/project/3015069.01/derived/BS-3/ses-mri03/gre/lowres/';
cd(base_folder)
out_base_folder = '/project/3015069.04/data/BS-3/';

fa_list = [5 10 20 30 40 50 70];

TE = linspace(2.12,54.31,18)*1e-3;
nb_TE = length(TE);

B0 = 3;
voxelSize = [1.4 1.4 1.4];

for kOri = 1:10
    orient_folder = [base_folder 'orientation-' num2str(kOri) '/'];
    out_orient_folder = [out_base_folder 'orientation-' num2str(kOri) '/'];
    mkdir(out_orient_folder);
    
    mask = load_nii_img_only([orient_folder 'BS-3_ses-mri03_acq-lowres_orientation-' num2str(kOri) '_gre_gradunwarp_mask_sphere.nii.gz']);
    matrixSize = size(mask);

    mask_replicate = repmat(mask, [1 1 1 nb_TE]);
    for fa = fa_list
        fa_folder = [orient_folder 'fa-' num2str(fa) '/'];
        out_fa_folder = [out_orient_folder 'fa-' num2str(fa) '/'];
        mkdir(out_fa_folder);
        
        magn_nii  = load_untouch_nii([fa_folder 'BS-3_ses-mri03_acq-lowres_orientation-' num2str(kOri) '_fa-' num2str(fa) '_part-mag_gre_gradunwarp.nii.gz']);
        magn = magn_nii.img;
        
        phase_nii  = load_untouch_nii([fa_folder 'BS-3_ses-mri03_acq-lowres_orientation-' num2str(kOri) '_fa-' num2str(fa) '_part-phase_gre_gradunwarp.nii.gz']);
        phase = phase_nii.img;
                
        %% phase unwrapping
        [~,~,unwarped_phase] = estimateTotalField(phase, magn, matrixSize,voxelSize,...
            'Method','Optimum weights','Unwrap','bestpath3d',...
            'unit','Hz','TE',TE,'B0',B0,...
            'mask',mask);
        
        phase_nii.img = single(unwarped_phase).*single(mask_replicate);
        save_untouch_nii(phase_nii, [out_fa_folder 'BS-3_ses-mri03_acq-lowres_orientation-' num2str(kOri) '_fa-' num2str(fa) '_part-phase_gre_gradunwarp_unwarped_masked.nii.gz']);
 
        phasemagn_nii_nii.img = single(unwarped_phase).*single(mask_replicate);
        save_untouch_nii(phase_nii, [out_fa_folder 'BS-3_ses-mri03_acq-lowres_orientation-' num2str(kOri) '_fa-' num2str(fa) '_part-phase_gre_gradunwarp_unwarped_masked.nii.gz']);
        
    end
end
