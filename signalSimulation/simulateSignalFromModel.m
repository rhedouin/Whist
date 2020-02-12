function [signal, field] = simulateSignalFromModel(axon_collection, model_parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%% Inputs required
% %%%% axon_collection is a structure of the white matter model where each
% element represent an axon with one required field
% % - data which corresponds to the myelin sheath
% % - optional fields: Centroid ,gRatio, axonEquiDiameter, myelinThickness
%
% %%%% model_parameters is a structure containing the parameters of the 3
% compartments as well as the simulation scan parameter
% % - intra_axonal : T2, weight 
% % - extra_axonal : T2, weight
% % - myelin : T2, weight, xi, xa
% The weight of each compartment reflects the water signal received, it comports 
% the water proton density and the T1 effect. 
% If you ignore the T1 effect the weight should be set to the proton density. 
% If you consider the T1 effect, you can call the following function to
% compute the weights, computeCompartmentSignalWeight(model_parameters).
%
% xi, xa are isotropic and anisotropic susceptibilities of myelin compare
% to intra_axonal and extra_axonal compartment which are set to 0
%
% % B0, MRI field strength (in Tesla)
% % field_direction, main magnetic field orientation relative to the model
% % TE, list of echo times (in second)
% % mask, area where the signal is estimated
%
% %%%%%%%%%%%%%% Outputs 
% %%%% signal is a structure containing multi GRE signal from each
% compartment and the total
% %%%% field is the field perturbation simulated from the WM model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


gamma = 42.6;
dims = size(model_parameters.mask);
number_dims = ndims(model_parameters.mask);

if  number_dims == 2
    [tensor_X, model]  = create2DTensorXFromAxonList(axon_collection, dims, model_parameters.myelin.xa, model_parameters.myelin.xi);
    
    field_complex = createFieldFrom2DTensorX(tensor_X, model_parameters.B0, gamma, model_parameters.field_direction);
    field = real(field_complex);
elseif number_dims == 3
    [tensor_X, model]  = create3DTensorXFromAxonList(axon_collection, dims, model_parameters.myelin.xa, model_parameters.myelin.xi);

    field_complex = createFieldFrom3DTensorX(tensor_X, model_parameters.B0, gamma, model_parameters.field_direction);
    field = real(field_complex);
    
else
    error('unexpected number of dimension, should be 2 or 3');
end

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
keyboard;
signal.intra_axonal = model_parameters.intra_axonal.weight * exp(-model_parameters.TE / model_parameters.intra_axonal.T2).* signal.intra_axonal / nb_pixel;
signal.myelin = model_parameters.myelin.weight * exp(-model_parameters.TE / model_parameters.myelin.T2).* signal.myelin / nb_pixel;
signal.extra_axonal = model_parameters.extra_axonal.weight * exp(-model_parameters.TE / model_parameters.extra_axonal.T2).*signal.extra_axonal / nb_pixel;

signal.total = signal.intra_axonal + signal.myelin + signal.extra_axonal;

end










