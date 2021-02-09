function corrected_field = tensorLorentzianCavityCorrection(tensor, model, model_parameters)
% Apply LorentzianCavityCorrection To Field Perturbation (see Biophysical mechanisms of phase contrast in gradient echo MRI, He 2009)

gamma = 42.6;
myelin_field_correction = model_parameters.myelin.xa * ((model_parameters.theta)^2 - 1/3);

corrected_tensor = tensor - myelin_field_correction * gamma * model_parameters.B0 .* (model == 1);
keyboard;
end