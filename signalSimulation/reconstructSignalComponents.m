function signal = reconstructSignalComponents(field, model, mask, time, current_FVF, current_gRatio, xi, xa, current_dir, dispersion)

% Compute signal components
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
    
    signal.Axon(l) = sum(C_axon, 'all');
    signal.Myelin(l) = sum(C_myelin, 'all');
    signal.Extra(l) = sum(C_extra, 'all');
end

signal.Axon = signal.Axon / nb_pixel;
signal.Myelin = signal.Myelin / nb_pixel;
signal.Extra = signal.Extra / nb_pixel;

signal.FVF = current_FVF;
signal.gRatio  = current_gRatio;
signal.xi = xi;
signal.xa = xa;
signal.current_dir = current_dir;
signal.dispersion = dispersion;
signal.time = time;

end
