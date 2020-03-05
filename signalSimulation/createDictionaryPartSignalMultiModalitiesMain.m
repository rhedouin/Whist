
clear
% close all

req_mem   = 4e9;
req_etime = 4000;

base_folder = '/project/3015069.04/';
job_folder = [base_folder 'temp/Jobs'];
signal_folder = [base_folder 'signal_components/multi_orientations/Porcine2/lowres/fix_xa_large_FVF_20_directions_tensor_mask/'];
dico_folder = [base_folder 'dictionaries/multi_orientations/Porcine2/lowres/fix_xa_large_FVF_20_directions_tensor_mask/'];

cd(job_folder)

FVFRange = (10 : 10 : 80);

T2MyelinRange = [4 8 12 16 20]*1e-3;
T2IntraExtraAxonalRange = (20 : 20 : 100)* 1e-3;
T2ExtraAxonalRange = (20 : 20 : 100)* 1e-3;

weightRange = [0.5 1 1.5 2 2.5 3];

nb_TE = 18;
nb_orientations = 6;

noise = 0;
nb_replic = 8;

it = 0;
experience_name = 'Porcine2';

options.include_theta = 1;

options.coordinate.classic_polar = 0;
options.coordinate.classic_cartesian = 0;
options.coordinate.polyfit_polar = 0;
options.coordinate.polyfit_cartesian = 1;
options.coordinate.polyfit_cartesian_demean = 0;

replic_list = [1, 2, 3, 4, 5, 7, 8];
for FVF = FVFRange
    for num = 1 : nb_replic
        it = it + 1;
        suffix = ['_train' num2str(num)];
        output_folder = [dico_folder '/dictionary_part/FVF' num2str(FVF) '_N400_train' num2str(num)];
        mkdir(output_folder)
        
        FVF_folder = [signal_folder 'FVF' num2str(FVF) '_N400' suffix '/'];
        signal_path = [FVF_folder 'Signal_FVF' num2str(FVF) suffix '.mat'];
        tic()

%         createDictionaryPartSignalMultiModalities(signal_path, output_folder, experience_name, T2MyelinRange, T2IntraExtraAxonalRange, weightRange, nb_TE, noise, FVF, nb_orientations, num, options);
        job{it} = qsubfeval(@createDictionaryPartSignalMultiModalities, signal_path, output_folder, experience_name, T2MyelinRange, T2IntraExtraAxonalRange, weightRange, nb_TE, noise, FVF, nb_orientations, num, options, 'memreq',  req_mem,  'timreq',  req_etime);
    end
end






