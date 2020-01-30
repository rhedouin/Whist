function [signal, field] = simulateSignalFromModel(axon_collection, mask, xa, xi,  T2out, T2myel, weight, time, B0, field_direction)

gamma = 42.6;
dims = size(mask);

[tensor_X, model]  = create2DTensorXFromAxonList(axon_collection, dims, xa, xi);

field_complex = createFieldFrom2DTensorX(tensor_X, B0, gamma, field_direction);
field = real(field_complex);

N = length(time);

if ndims(field) == 3
    model = repmat(model, 1, 1, size(field,3));
    mask = repmat(mask, 1, 1, size(field,3));
    display('3D field')
end

model(mask == 0) = -1;

omega = 2*pi*field;

axon_index = (model== 0.5);
myelin_index = (model == 1);
extra_index = (model == 0 );

nb_pixel = sum(mask, 'all');
for l = 1:N
    
    C_axon = exp(1i*time(l)*omega(axon_index));
    C_myelin = exp(1i*time(l)*omega(myelin_index));
    C_extra = exp(1i*time(l)*omega(extra_index));
    
    signal_Axon(l) = sum(C_axon, 'all');
    signal_Myelin(l) = sum(C_myelin, 'all');
    signal_Extra(l) = sum(C_extra, 'all');
end

signal_Axon = exp(-time/T2myel).*signal_Axon / nb_pixel;
signal_Myelin = weight*exp(-time/T2out).*signal_Myelin / nb_pixel;
signal_Extra = exp(-time/T2out).*signal_Extra / nb_pixel;

signal = signal_Axon + signal_Myelin + signal_Extra;

end










