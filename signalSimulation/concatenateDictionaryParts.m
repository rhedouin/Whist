clear
% close all

base_folder = '/project/3015069.04/';
dico_folder = [base_folder 'dictionaries/multi_orientations/Porcine2/'];

FVFRange = (40 : 5 : 85);

lFVF = length(FVFRange);

nb_replica = 8;
nb_orientation = 6;
nb_TE = 12;
noise = '1';
experience_name = 'Porcine2';

suffix = 'polyfit_cartesian_without_theta';

tic()
display('Concatenation ...')

input_folder = [dico_folder 'dictionary_part/FVF40_N400_train1/'];
dict_path = [input_folder 'SignalWithNoise'  noise '_FVF40_replic1_' num2str(nb_orientation) 'orientations_' num2str(nb_TE) 'TE_' experience_name '_fix_xa_' suffix '.h5py'];

dims_signal = size(h5read(dict_path, '/SignalValues'));
news_dims_signal = [dims_signal(1) lFVF dims_signal(2:end) nb_replica];
news_dims_para = [lFVF dims_signal(2:end) nb_replica];

SignalValues = zeros(news_dims_signal);

gRatioValues = zeros(news_dims_para);
xiValues = zeros(news_dims_para);
xaValues = zeros(news_dims_para);
T2outValues = zeros(news_dims_para);
T2myelValues = zeros(news_dims_para);
weightValues = zeros(news_dims_para);
thetaValues = zeros(news_dims_para);
directionsValues = zeros([3, news_dims_para]);

for k = 1 : lFVF
    FVF = FVFRange(k)
    for num = 1 : nb_replica
        num
        input_folder = [dico_folder 'dictionary_part/FVF' num2str(FVF) '_N400_train' num2str(num) '/'];
        dict_path = [input_folder 'SignalWithNoise'  noise '_FVF' num2str(FVF) '_replic' num2str(num) '_' num2str(nb_orientation) 'orientations_' num2str(nb_TE) 'TE_Porcine2_fix_xa_' suffix '.h5py'];
        
        SignalValues(:, k, :, :, :, :, :, :, :,  num) = h5read(dict_path, '/SignalValues');
        
        gRatioValues(k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/gRatioValues');
        xiValues(k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/xiValues');
        xaValues(k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/xaValues');
        T2outValues(k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/T2outValues');
        T2myelValues(k, :, :, :, :, :, :, :, num) = h5read(dict_path, '/T2myelValues');
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
xiRange = h5read(dict_path, '/xiRange');
xaRange = h5read(dict_path, '/xaRange');
T2outRange = h5read(dict_path, '/T2outRange');
T2myelRange = h5read(dict_path, '/T2myelRange');
weightRange = h5read(dict_path, '/weightRange');

time = h5read(dict_path, '/time');
sphere_rotations = h5read(dict_path, '/sphere_rotations');

infoDico = 'In order, FVF, gRatio, xi, dir, T2myel, T2out, weight, replic';
infoSignal = ['concatenation of ' num2str(nb_orientation) ' rotations each composed of the theta angle (B0 angle) 12 TE signal real and imaginary polyfit normalized'];

prefix_name = ['SignalWithNoise' noise];

base_name = [prefix_name '_' num2str(nb_replica) 'rep_' num2str(nb_orientation) 'orientations_' num2str(nb_TE) 'TE_' experience_name '_fix_xa_' suffix];
signal_name = [base_name '.h5py'];

display('Save ...')
tic()
save([dico_folder signal_name], 'SignalValues', ...
    'FVFRange', 'gRatioRange', 'xiRange', 'xaRange', 'T2myelRange','T2outRange','weightRange', 'FVFValues', ...
    'gRatioValues', 'xiValues', 'xaValues', 'T2myelValues', 'T2outValues', 'directionsValues', 'thetaValues', ...
    'weightValues', 'nb_replica', 'time', 'infoDico', 'infoSignal', 'sphere_rotations', '-v7.3')
toc()
display('Done')












