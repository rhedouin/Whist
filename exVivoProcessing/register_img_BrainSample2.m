clear;

base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/';
parameter_folder = [base_folder 'concatenate_signals_9_orientations/parameter_maps/BrainSample2LorentzinaCorrection/'];
ref_folder = [base_folder 'references/'];

input = [ref_folder 'BrainSample-2_ses-03_mp2rage_orientation-4_t1_brain_2_BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1_brain'];
ref = [ref_folder 'ref'];
% other = {'BrainSample-2_ses-03_FVF_mean_20_directions_polyfit_cartesian_with_theta_mask_noise4'};
ext = '_register.nii.gz';

fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60', 'mean', 'std'}
% fa_list = {'mean', 'std'}

parameter_list = {'FVF', 'gRatio', 'T2Myelin', 'T2IntraExtraAxonal', 'weight', 'xiMyelin'};
% parameter_list = { 'R2Myelin', 'R2IntraExtraAxonal'};

other = {};

mask = load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_gre_orientation-4_brain_mask.nii.gz');

for k = 1:length(fa_list)
    fa = fa_list{k}
    fa_folder = [parameter_folder fa '/'];
    
    for l = 1:length(parameter_list)
        parameter = parameter_list{l};
        other{end+1} = [fa_folder 'BrainSample-2_ses-03_' parameter '_' fa '_polyfit_cartesian_with_theta_noise4'];
    end
end
    
  
% other{end+1} = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/references/BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1';

flirtAB2C(input, other, ref, ext)