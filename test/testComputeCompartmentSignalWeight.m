% model_parameters = computeCompartmentSignalWeight(model_parameters);
% clear

model_parameters.myelin.T2 = 15*1e-3;
model_parameters.myelin.T1 = 500*1e-3;
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

model_parameters.flip_angle = 20;
model_parameters.TR = 60*1e-3;
% 
model_parameters.include_proton_density = 1;
model_parameters.include_T1_effect = 1;

list_fa = [5, 10, 25, 40, 70];


for l = 1:length(list_fa)
    fa = list_fa(l);
    
    model_parameters.flip_angle = fa;
    
    model_parameters = computeCompartmentSignalWeight(model_parameters);
    relative_weight( = model_parameters.myelin.weight / model_parameters.intra_axonal.weight;
end