function signal_components = reconstructSignalComponents(field, model, model_parameters)

% Compute signal components
N = length(model_parameters.TE);

if ndims(field) == 3
    model = repmat(model, 1, 1, size(field,3));
    model_parameters.mask = repmat(model_parameters.mask, 1, 1, size(field,3));
    display('3D field')
end

model(model_parameters.mask == 0) = -1;

omega = 2*pi*field;

axon_index = (model== 0.5);
myelin_index = (model == 1);
extra_index = (model == 0 );

nb_pixel = sum(model_parameters.mask, 'all');

for l = 1:N
    
    C_axon = exp(1i*model_parameters.TE(l)*omega(axon_index));
    C_myelin = exp(1i*model_parameters.TE(l)*omega(myelin_index));
    C_extra = exp(1i*model_parameters.TE(l)*omega(extra_index));
    
    signal_components.intra_axonal.signal(l) = sum(C_axon, 'all');
    signal_components.myelin.signal(l) = sum(C_myelin, 'all');
    signal_components.extra_axonal.signal(l) = sum(C_extra, 'all');
end

signal_components.intra_axonal.signal = signal_components.intra_axonal.signal / nb_pixel;
signal_components.myelin.signal = signal_components.myelin.signal / nb_pixel;
signal_components.extra_axonal.signal = signal_components.extra_axonal.signal / nb_pixel;

signal_components.FVF = model_parameters.FVF;
signal_components.g_ratio  = model_parameters.g_ratio;
signal_components.myelin.xi = model_parameters.myelin.xi;
signal_components.myelin.xa = model_parameters.myelin.xa;
signal_components.current_dir = model_parameters.current_dir;
signal_components.dispersion = model_parameters.dispersion;
signal_components.TE = model_parameters.TE;

if isfield(model_parameters, 'intra_axonal')
    signal_components.intra_axonal.xi = model_parameters.intra_axonal.xi;
else
    signal_components.intra_axonal.xi = 0;
end
if isfield(model_parameters, 'extra_axonal')
    signal_components.extra_axonal.xi = model_parameters.extra_axonal.xi;
else
    signal_components.extra_axonal.xi = 0;
end

end
