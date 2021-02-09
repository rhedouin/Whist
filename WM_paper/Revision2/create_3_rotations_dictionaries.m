clear
dict_folder = '/project/3015069.04/dictionaries/multi_orientations/BrainSample2LorentzianCorrection/';
cd(dict_folder)
input_dict_name = 'SignalWithNoise4_8rep_9rotations_12TE_BrainSample2LorentzianCorrection_polyfit_cartesian_with_theta.h5py';

% h5disp(input_dict_name);
FVFRange = h5read(input_dict_name, '/FVFRange');
FVFValues = h5read(input_dict_name, '/FVFValues');
T2IntraExtraAxonalRange = h5read(input_dict_name, '/T2IntraExtraAxonalRange');
T2IntraExtraAxonalValues = h5read(input_dict_name, '/T2IntraExtraAxonalValues');
T2MyelinRange = h5read(input_dict_name, '/T2MyelinRange');
T2MyelinValues = h5read(input_dict_name, '/T2MyelinValues');
TE = h5read(input_dict_name, '/TE');
directionsValues = h5read(input_dict_name, '/directionsValues');
gRatioRange = h5read(input_dict_name, '/gRatioRange');
gRatioValues = h5read(input_dict_name, '/gRatioValues');
infoDico = h5read(input_dict_name, '/infoDico');
infoSignal = h5read(input_dict_name, '/infoSignal');
nb_replica = h5read(input_dict_name, '/nb_replica');
sphere_rotations = h5read(input_dict_name, '/sphere_rotations');
thetaValues = h5read(input_dict_name, '/thetaValues');
weightRange = h5read(input_dict_name, '/weightRange');
weightValues = h5read(input_dict_name, '/weightValues');
xaMyelinRange = h5read(input_dict_name, '/xaMyelinRange');
xaMyelinValues = h5read(input_dict_name, '/xaMyelinValues');
xiMyelinRange = h5read(input_dict_name, '/xiMyelinRange');
xiMyelinValues = h5read(input_dict_name, '/xiMyelinValues');

SignalValues_save = h5read(input_dict_name, '/SignalValues');

list_orientations = {[4 5 7], [1 5 7], [2 4 6], [2 4 5], [6 7 9], [2 6 9], ...
                     [4 6 7], [5 7 9], [1 2 5], [2 3 4]};

orient_unit = 25;
for k = 1:length(list_orientations)
    current_orient = list_orientations{k};
    subset_orient = [(current_orient(1) - 1) *orient_unit + 1 : current_orient(1) * orient_unit, ...
                     (current_orient(2) - 1) *orient_unit + 1 : current_orient(2) * orient_unit, ...
                     (current_orient(3) - 1) *orient_unit + 1 : current_orient(3) * orient_unit]

    SignalValues = SignalValues_save(subset_orient, :, :, :, :, :, :, :, :, :);
    
    output_dict_name = ['SignalWithNoise4_8rep_' num2str(current_orient(1)) num2str(current_orient(2)) num2str(current_orient(3)) 'rotations_12TE_BrainSample2LorentzianCorrection_polyfit_cartesian_with_theta.h5py'];
    save([dict_folder output_dict_name], 'SignalValues', ...
        'FVFRange', 'gRatioRange', 'xiMyelinRange', 'xaMyelinRange', 'T2MyelinRange','T2IntraExtraAxonalRange','weightRange', 'FVFValues', ...
        'gRatioValues', 'xiMyelinValues', 'xaMyelinValues', 'T2MyelinValues', 'T2IntraExtraAxonalValues', 'directionsValues', 'thetaValues', ...
        'weightValues', 'nb_replica', 'TE', 'infoDico', 'infoSignal', 'sphere_rotations', '-v7.3');

end





