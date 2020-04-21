% Separate concatenate signal

base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/';
concatenate_folder = [base_folder 'concatenate_signals/'];

fa_list = {'05', '10', '15', '20', '35', '60'};
lRot = 25;

nb_rotations = 9;

for k = 1:length(fa_list)
    fa = fa_list{k};
    
    signal_path = [concatenate_folder 'BrainSample-2_ses-03_all_orientations_fa-' fa '_concatenate_signal_polyfit_cartesian_with_theta_2_orientation-4.nii.gz'];
    signal_nii = load_untouch_nii(signal_path);
    
    signal_nii.hdr.dime.dim(5) = lRot;
    signal = signal_nii.img;
    
    for rot=1:nb_rotations
        
        orientation_folder = [base_folder 'orientation' num2str(rot) '/'];
        mkdir(orientation_folder);
        
        signal_one_rotation = signal(:, :, :, (rot-1) *lRot + 1: rot*lRot);        
        signal_nii.img = signal_one_rotation;
        
        save_untouch_nii(signal_nii, [orientation_folder 'BrainSample-2_ses-03_orientation' num2str(rot) '_fa-' fa '_concatenate_signal_polyfit_cartesian_with_theta_2_orientation-4.nii.gz']);
    end
end


