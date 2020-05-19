list_orientations = [157, 246, 245, 679, 269, 467, 579, 125, 234];
list_fa = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60'};

base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/';

lVec = 25;
for k = 1:length(list_fa)
    fa = list_fa{k}
    signal_path = [base_folder 'concatenate_signals_9_orientations/BrainSample-2_ses-03_all_orientations_' fa '_concatenate_signal_polyfit_cartesian_with_theta_2_orientation-4.nii.gz'];
    
    for l = 1:length(list_orientations)
        l
        input_nii = load_untouch_nii(signal_path);
        output_nii = input_nii;
        current_orientation = list_orientations(l);
        current_folder = [base_folder 'concatenate_signals_3_orientations/' num2str(current_orientation) 'orientations/'];
        mkdir(current_folder);
        
        num(3) = mod(current_orientation, 10);
        num(2) = (mod(current_orientation, 100) - num(3))/10;
        num(1) = (current_orientation - num(2)*10 - num(3))/100;
        
        list_signal = [];
        for m = 1:3      
            list_signal = [list_signal (num(m) - 1) * lVec + 1 : num(m) * lVec];
        end
        
        output_nii.img = input_nii.img(:, :, :, list_signal);
        output_nii.hdr.dime.dim(5) = 3*lVec;
        
        save_untouch_nii(output_nii, [current_folder 'BrainSample-2_ses-03_' num2str(current_orientation) '_orientations_' fa '_concatenate_signal_polyfit_cartesian_with_theta_2_orientation-4.nii.gz']);
    end
end


