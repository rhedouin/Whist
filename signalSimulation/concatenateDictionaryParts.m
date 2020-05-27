function out = concatenateDictionaryParts(dict_folder, FVFRange, prefix_name, nb_replica, nb_rotations, nb_TE, experience_name, suffix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function concatenate dictionary parts obtained from
% computeTotalSignalDictionaryPart. The final dictionary obtained is saved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic()
display('Concatenation ...')

lFVF = length(FVFRange);

replic_list = 1:nb_replica;
temp_FVF = FVFRange(1);
input_folder = [dict_folder 'dictionary_part/FVF' num2str(temp_FVF) '_N400_train1/'];

dict_path = [input_folder prefix_name '_FVF' num2str(temp_FVF) '_replic1_' num2str(nb_rotations) 'rotations_' num2str(nb_TE) 'TE_' experience_name suffix '.h5py'];

dims_signal = size(h5read(dict_path, '/SignalValues'));
news_dims_signal = [dims_signal(1) lFVF dims_signal(2:end) nb_replica];
news_dims_para = [lFVF dims_signal(2:end) nb_replica];

SignalValues = zeros(news_dims_signal, 'single');

gRatioValues = zeros(news_dims_para, 'single');
xiMyelinValues = zeros(news_dims_para, 'single');
xaMyelinValues = zeros(news_dims_para, 'single');
T2IntraExtraAxonalValues = zeros(news_dims_para, 'single');
T2MyelinValues = zeros(news_dims_para, 'single');
weightValues = zeros(news_dims_para, 'single');
thetaValues = zeros(news_dims_para, 'single');
directionsValues = zeros([3, news_dims_para], 'single');


for k = 1 : lFVF
    FVF = FVFRange(k)
    for l = 1 : length(replic_list)
        num = replic_list(l)
        input_folder = [dict_folder 'dictionary_part/FVF' num2str(FVF) '_N400_train' num2str(num) '/'];
        dict_path = [input_folder prefix_name '_FVF' num2str(FVF) '_replic' num2str(num) '_' num2str(nb_rotations) 'rotations_' num2str(nb_TE) 'TE_' experience_name suffix '.h5py'];
        
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
infoSignal = ['concatenation of ' num2str(nb_rotations) ' rotations each composed of the theta angle (B0 angle) 12 TE signal real and imaginary polyfit normalized'];

base_name = [prefix_name '_' num2str(nb_replica) 'rep_' num2str(nb_rotations) 'rotations_' num2str(nb_TE) 'TE_' experience_name suffix];
signal_name = [base_name '.h5py'];


display('Save ...')
tic()
save([dict_folder signal_name], 'SignalValues', ...
    'FVFRange', 'gRatioRange', 'xiMyelinRange', 'xaMyelinRange', 'T2MyelinRange','T2IntraExtraAxonalRange','weightRange', 'FVFValues', ...
    'gRatioValues', 'xiMyelinValues', 'xaMyelinValues', 'T2MyelinValues', 'T2IntraExtraAxonalValues', 'directionsValues', 'thetaValues', ...
    'weightValues', 'nb_replica', 'TE', 'infoDico', 'infoSignal', 'sphere_rotations', '-v7.3')
toc()
display('Done')

end









