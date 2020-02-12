function model_parameters = computeCompartmentSignalWeight(model_parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%% Inputs required
% model parameters with 
% T1 and the proton density of each compartment
% TR
% flip angle 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E1 = exp(-model_parameters.TR/model_parameters.myelin.T1);
model_parameters.myelin.weight = model_parameters.myelin.proton_density * sind(model_parameters.flip_angle) * (1 - E1) / (1 - cosd(model_parameters.flip_angle) * E1);

E1 = exp(-model_parameters.TR/model_parameters.intra_axonal.T1);
model_parameters.intra_axonal.weight = model_parameters.intra_axonal.proton_density * sind(model_parameters.flip_angle) * (1 - E1) / (1 - cosd(model_parameters.flip_angle) * E1);

E1 = exp(-model_parameters.TR/model_parameters.extra_axonal.T1);
model_parameters.extra_axonal.weight = model_parameters.extra_axonal.proton_density * sind(model_parameters.flip_angle) * (1 - E1) / (1 - cosd(model_parameters.flip_angle) * E1);

end