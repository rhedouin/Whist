function [signal, field] = simulateSignalFromModel(axon_collection, model_parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% axon_collection 



gamma = 42.6;
dims = size(model_parameters.mask);

[tensor_X, model]  = create2DTensorXFromAxonList(axon_collection, dims, model_parameters.myelin.xa, model_parameters.myelin.xi);

field_complex = createFieldFrom2DTensorX(tensor_X, model_parameters.B0, gamma, model_parameters.field_direction);
field = real(field_complex);

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

end










