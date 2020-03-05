clear
% close all

base_folder = '/project/3015069.04/';
dico_folder = [base_folder 'dictionaries/multi_orientations/Porcine2/lowres/fix_xa_large_FVF_20_directions_tensor_mask/'];

FVFRange = (10 : 10 : 80);

lFVF = length(FVFRange);

nb_replica = 8;
nb_orientation = 6;
nb_TE = 18;
noise = '0';
experience_name = 'Porcine2';

input_suffix = 'fix_xa_polyfit_cartesian_with_theta';
output_suffix = 'fix_xa_large_FVF_20_directions_polyfit_cartesian_with_theta';

tic()
display('Concatenation ...')

input_folder = [dico_folder 'dictionary_part/FVF40_N400_train1/'];
dict_path = [input_folder 'SignalWithNoise'  noise '_FVF40_replic1_' num2str(nb_orientation) 'orientations_' num2str(nb_TE) 'TE_' experience_name '_' input_suffix '.h5py'];

dims_signal = size(h5read(dict_path, '/SignalValues'));
news_dims_signal = [dims_signal(1) lFVF dims_signal(2:end) nb_replica];
news_dims_para = [lFVF dims_signal(2:end) nb_replica];

SignalValues = zeros(news_dims_signal);

gRatioValues = zeros(news_dims_para);
xiMyelinValues = zeros(news_dims_para);
xaMyelinValues = zeros(news_dims_para);
T2IntraExtraAxonalValues = zeros(news_dims_para);
T2MyelinValues = zeros(news_dims_para);
weightValues = zeros(news_dims_para);
thetaValues = zeros(news_dims_para);
directionsValues = zeros([3, news_dims_para]);

for k = 1 : lFVF
    FVF = FVFRange(k)
    for num = 1 : nb_replica
        num
        input_folder = [dico_folder 'dictionary_part/FVF' num2str(FVF) '_N400_train' num2str(num) '/'];
        dict_path = [input_folder 'SignalWithNoise'  noise '_FVF' num2str(FVF) '_replic' num2str(num) '_' num2str(nb_orientation) 'orientations_' num2str(nb_TE) 'TE_' experience_name '_' input_suffix '.h5py'];
        
        SignalValues(:, k, :, :, :, :, :, :, :,  num) = h5read(dict_path, '/SignalValues');
        
        gRatioValues(k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/gRatioValues');
        xiMyelinValues(k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/xiMyelinValues');
        xaMyelinValues(k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/xaMyelinValues');
        T2IntraExtraAxonalValues(k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/T2IntraExtraAxonalValues');
        T2MyelinValues(k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/T2MyelinValues');
        weightValues(k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/weightValues');
        thetaValues(k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/thetaValues');
        directionsValues(:, k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/directionsValues');      
        
    end
end

toc()
display('Done ...')
dims = size(gRatioValues);
FVFValues = repmat(FVFRange', [1 dims(2:end)])./100;

gRatioRange = h5read(dict_path, '/gRatioRange');
xiMyelinRange = h5read(dict_path, '/xiMyelinRange');
xaMyelinRange = h5read(dict_path, '/xaMyelinRange');
T2IntraExtraAxonalRange = h5read(dict_path, '/T2IntraExtraAxonalRange');
T2MyelinRange = h5read(dict_path, '/T2MyelinRange');
weightRange = h5read(dict_path, '/weightRange');

TE = h5read(dict_path, '/TE');
sphere_rotations = h5read(dict_path, '/sphere_rotations');

infoDico = 'In order, FVF, gRatio, xiMyelin, xaMyelin, dir, T2Myelin, T2IntraExtraAxonalRange, weight, replic';
infoSignal = ['concatenation of ' num2str(nb_orientation) ' rotations each composed of the theta angle (B0 angle) 12 TE signal real and imaginary polyfit normalized'];

prefix_name = ['SignalWithNoise' noise];

base_name = [prefix_name '_' num2str(nb_replica) 'rep_' num2str(nb_orientation) 'orientations_' num2str(nb_TE) 'TE_' experience_name '_' output_suffix '.h5py'];
signal_name = [base_name '.h5py'];

display('Save ...')
tic()
save([dico_folder signal_name], 'SignalValues', ...
    'FVFRange', 'gRatioRange', 'xiMyelinRange', 'xaMyelinRange', 'T2MyelinRange','T2IntraExtraAxonalRange','weightRange', 'FVFValues', ...
    'gRatioValues', 'xiMyelinValues', 'xaMyelinValues', 'T2MyelinValues', 'T2IntraExtraAxonalValues', 'directionsValues', 'thetaValues', ...
    'weightValues', 'nb_replica', 'TE', 'infoDico', 'infoSignal', 'sphere_rotations', '-v7.3')
toc()
display('Done')












