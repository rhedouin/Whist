cd '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/test/output'

A = load_untouch_nii('Sepia_total-field.nii.gz')
total_field = A.img;

A = load_untouch_nii('Sepia_local-field.nii.gz')
local_field = A.img;

A = load_untouch_nii('Sepia_unwrapped-phase.nii.gz')
phase = A.img;

load('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/test/BrainSample-2_ses-03_gre_orientation-1_fa-05_header.mat')

back_field = total_field - local_field;
back_field_rep = repmat(back_field, [1 1 1 12]);

dims = size(local_field);
TE_rep = permute(repmat(TE(:), [1 dims(1) dims(2) dims(3)]), [2 3 4 1]);

phase_corrected_plus = phase .* (back_field_rep.*TE_rep);
phase_corrected_minus = phase .* (-back_field_rep.*TE_rep);

diff_minus = phase - phase_corrected_minus;

A.img = phase_corrected_minus;
save_untouch_nii(A, 'phase_corrected_minus.nii.gz')

A.img = phase_corrected_plus;
save_untouch_nii(A, 'phase_corrected_plus.nii.gz')

A.img = diff_minus;
save_untouch_nii(A, 'diff_minus.nii.gz')

A.hdr.dime.dim(1) = 3;
A.hdr.dime.dim(5) = 1;
A.img = back_field;
save_untouch_nii(A, 'back_field.nii.gz');
keyboard;

unwrap