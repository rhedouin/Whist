% test compute average weight
clear;
close all
data_folder = '/project/3015069.01/derived/Porcine-2/ses-mri01/concatenate_signals/lowres/individual_samples/parameter_maps/20_directions_polyfit_cartesian_with_theta/';

list_fa = [5, 10, 25, 40, 70];
nb_samples = 7;

figure
for k = 1:nb_samples
for l = 1:length(list_fa)

    fa = list_fa(l);
    weight_map = load_nii_img_only([data_folder 'Porcine-2_ses-mri01_small_sample_' num2str(k) '_fa-' num2str(fa) '_weight_20_directions_polyfit_cartesian_with_theta.nii.gz']);
    mask = load_nii_img_only(['/project/3015069.01/derived/Porcine-2/ses-mri01/concatenate_signals/lowres/individual_samples/Porcine-2_ses-mri01_small_sample_' num2str(k) '_mask.nii.gz']);
    
    weight(k,l) = sum(weight_map, 'all') /  sum(mask, 'all');
end
plot(list_fa, weight(k,:))
hold on
end


model_parameters.myelin.T2 = 15*1e-3;
model_parameters.myelin.T1 = 400*1e-3;
model_parameters.myelin.proton_density= 0.5; 

model_parameters.myelin.xi = -0.1;  % myelin anisotropic susceptibility (ppm)
model_parameters.myelin.xa = -0.1;  % myelin isotropic susceptibility (ppm)

% intra axonal (required: T2)
model_parameters.intra_axonal.T2 = 50*1e-3;
model_parameters.intra_axonal.T1 = 1.5;
model_parameters.intra_axonal.proton_density= 1; 
model_parameters.intra_axonal.xi= 0; 

% extra axonal (required: T2)
model_parameters.extra_axonal.T2 = 50*1e-3;
model_parameters.extra_axonal.T1 = 1.5;
model_parameters.extra_axonal.proton_density= 1; 
model_parameters.extra_axonal.xi= 0; 

model_parameters.TR = 60*1e-3;
% 
model_parameters.include_proton_density = 1;
model_parameters.include_T1_effect = 1;

list_fa = [5, 10, 25, 40, 70];

for l = 1:length(list_fa)
    fa = list_fa(l);
    
    model_parameters.flip_angle = fa;
    
    model_parameters = computeCompartmentSignalWeight(model_parameters);
    relative_weight(l) = model_parameters.myelin.weight / model_parameters.intra_axonal.weight;
end

plot(list_fa, relative_weight)


leg = legend('1', '2', '3', '4', '5', '6', '7', 'theoretical');
title(leg, 'sample')

