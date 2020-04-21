function [signal, field] = simulateSignalWithDispersionFromModel(axon_collection, model_parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%% Inputs required
% %%%% axon_collection is a structure of the white matter model where each
% element represent an axon with 
% % - required field: data which corresponds to the myelin sheath
% % - optional fields: Centroid ,gRatio, axonEquiDiameter, myelinThickness
%
% %%%% model_parameters is a structure containing the parameters of the 3
% compartments as well as the simulation scan parameter
% % - intra_axonal: T2, weight 
% % - extra_axonal: T2, weight
% % - myelin: T2, weight, xi, xa
%
% xi, xa are isotropic and anisotropic susceptibilities of myelin compare
% to a reference (intra/extra axonal compartment by defaut)
% weight reflects the water signal received, see computeCompartmentSignalWeight
%
% % B0, MRI field strength (in Tesla)
% % field_direction, main magnetic field orientation relative to the model
% % TE, list of echo times (in second)
% % mask, area where the signal is estimated
%
 %%%%%%%%%%%%%%% Inputs optional
% %%%% model_parameters
% % - intra_axonal, extra_axonal: xi (set to 0 by defauts)

% %%%%%%%%%%%%%% Outputs 
% %%%% signal is a structure containing the multi GRE signal from each
% compartment and the total signal (sum of compartments)
% %%%% field is the field perturbation simulated from the WM model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

number_dims = length(model_parameters.dims);

[tensor_X, model]  = create2DTensorXFromAxonList(axon_collection, model_parameters);

[orientations, dispersion] = generateVMFsample_for_dispersion(model_parameters.field_direction, model_parameters.kappa, model_parameters.dispersion, model_parameters.nb_orientation_for_dispersion);

for k = 1:length(orientations)
    model_parameters.field_direction = orientations(k,:);
    field_complex(:,:,k) = createFieldFrom2DTensorX(tensor_X, model_parameters);
end

field = real(field_complex);

model = repmat(model, 1, 1, size(field,3));
model_parameters.mask = repmat(model_parameters.mask, 1, 1, size(field,3));    

N = length(model_parameters.TE);
model(model_parameters.mask == 0) = -1;

omega = 2*pi*field;

intra_axonal_index = (model == 0.5);
myelin_index = (model == 1);
extra_axonal_index = (model == 0 );

nb_pixel = sum(model_parameters.mask, 'all');

for l = 1:N    
    C.intra_axonal = exp(1i*model_parameters.TE(l)*omega(intra_axonal_index));
    C.myelin = exp(1i*model_parameters.TE(l)*omega(myelin_index));
    C.extra_axonal = exp(1i*model_parameters.TE(l)*omega(extra_axonal_index));
    
    signal.intra_axonal(l) = sum(C.intra_axonal, 'all');
    signal.myelin(l) = sum(C.myelin, 'all');
    signal.extra_axonal(l) = sum(C.extra_axonal, 'all');
end

signal.intra_axonal = model_parameters.intra_axonal.weight * exp(-model_parameters.TE / model_parameters.intra_axonal.T2).* signal.intra_axonal / nb_pixel;
signal.myelin = model_parameters.myelin.weight * exp(-model_parameters.TE / model_parameters.myelin.T2).* signal.myelin / nb_pixel;
signal.extra_axonal = model_parameters.extra_axonal.weight * exp(-model_parameters.TE / model_parameters.extra_axonal.T2).*signal.extra_axonal / nb_pixel;

signal.total = signal.intra_axonal + signal.myelin + signal.extra_axonal;
signal.total_normalized = signal.total  / abs(signal.total(1));
end










