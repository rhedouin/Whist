function corrected_field = applyFieldLorentzianCavityCorrection(field, susceptibility_Z, model_parameters)
% Apply LorentzianCavityCorrection To Field Perturbation (see Biophysical mechanisms of phase contrast in gradient echo MRI, He 2009)

gamma = 42.6; % MHz

myelin_field_correction = susceptibility_Z * (cos(model_parameters.theta)^2 - 1/3) / 2;
corrected_field = field - gamma * model_parameters.B0 * myelin_field_correction;

end