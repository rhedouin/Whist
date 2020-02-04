% subsample dictionary orientations
nb_orientation = 15;

base_folder = '/project/3015069.04/dictionaries/multi_orientations/theorically_good_15_rotations/';
input_path = [base_folder 'SignalWithNoise05_5rep_' num2str(nb_orientation) 'orientations_fix_xa_polyfit_cartesian_with_theta.h5py'];

FVFRange = h5read(input_path, '/FVFRange');
FVFValues = h5read(input_path, '/FVFValues');
T2myelRange = h5read(input_path, '/T2myelRange');
T2myelValues = h5read(input_path, '/T2myelValues');
T2outRange = h5read(input_path, '/T2outRange');
T2outValues = h5read(input_path, '/T2outValues');
directionsValues = h5read(input_path, '/directionsValues');
gRatioRange = h5read(input_path, '/gRatioRange');
gRatioValues = h5read(input_path, '/gRatioValues');
infoDico = h5read(input_path, '/infoDico');
infoSignal = h5read(input_path, '/infoSignal');
nb_replica = h5read(input_path, '/nb_replica');
sphere_rotations = h5read(input_path, '/sphere_rotations');
thetaValues = h5read(input_path, '/thetaValues');
time = h5read(input_path, '/time');
weightRange = h5read(input_path, '/weightRange');
weightValues = h5read(input_path, '/weightValues');
xaRange = h5read(input_path, '/xaRange');
xaValues = h5read(input_path, '/xaValues');
xiRange = h5read(input_path, '/xiRange');
xiValues = h5read(input_path, '/xiValues');

display('Load signal ...')
SignalValues = h5read(input_path, '/SignalValues');
display('Done')

for k = 1:nb_orientation - 1
    display(['nb orientation : ' num2str(k)])

    infoSignal = ['concatenation of ' num2str(k) ' rotations each composed of the theta angle (B0 angle) 12 TE signal real and imaginary polyfit normalized'];
   
    SignalValues = SignalValues(1:end-25, :, :, :, :, :, :, :, :);

    output_path = [base_folder 'SignalWithNoise05_5rep_' num2str(nb_orientation - k) 'orientations_fix_xa_polyfit_cartesian_with_theta.h5py'];
    
    save(output_path, 'SignalValues', ...
    'FVFRange', 'gRatioRange', 'xiRange', 'xaRange', 'T2myelRange','T2outRange','weightRange', 'FVFValues', ...
    'gRatioValues', 'xiValues', 'xaValues', 'T2myelValues', 'T2outValues', 'directionsValues', 'thetaValues', ...
    'weightValues', 'nb_replica', 'time', 'infoDico', 'infoSignal', 'sphere_rotations', '-v7.3')    
end
