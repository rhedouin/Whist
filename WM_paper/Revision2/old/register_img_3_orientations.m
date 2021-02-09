clear;

base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/';

input = [base_folder 'references/BrainSample-2_ses-03_mp2rage_orientation-4_t1_brain_2_BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1_brain'];
ref = [base_folder 'references/ref'];
signal_folder = [base_folder 'concatenate_signals_3_orientations/'];
ext = '_register.nii.gz';

fa_list = {'fa-05', 'fa-10', 'fa-15', 'fa-20', 'fa-35', 'fa-60', 'mean', 'std'};
parameter_list = {'FVF', 'gRatio', 'T2Myelin', 'T2IntraExtraAxonal', 'weight', 'xiMyelin'};
list_orientations = {'457', '157', '246', '245', '679', '269', '467', '579', '125', '234'};
    
other = {};

mask = load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_gre_orientation-4_brain_mask.nii.gz');

for m = 1 : length(list_orientations)
    current_orientation = list_orientations{m};
    parameter_folder = [signal_folder current_orientation 'orientations/parameter_maps/BrainSample2LorentzinaCorrection/'];
    for k = 1:length(fa_list)
        fa = fa_list{k};
        
        for l = 1:length(parameter_list)
            parameter = parameter_list{l};
            other{end+1} = [parameter_folder 'BrainSample-2_ses-03_' parameter '_' fa '_polyfit_cartesian_with_theta_noise4'];
        end
    end
end
    
  
other{end+1} = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/references/BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1';

flirtAB2C(input, other, ref, ext)