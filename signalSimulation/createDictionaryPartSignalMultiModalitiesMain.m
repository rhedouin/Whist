clear
% close all

req_mem   = 5e9;
req_etime = 5000;

base_folder = '/project/3015069.04/';
job_folder = [base_folder 'temp/Jobs'];
signal_folder = [base_folder 'signal_components/multi_orientations/Porcine2/without_dispersion/'];
dico_folder = [base_folder 'dictionaries/multi_orientations/Porcine2/'];

cd(job_folder)

FVFRange = (40 : 5 : 85);

T2myelRange = [4 8 12 16 20]*1e-3;
T2outRange = (20 : 20 : 100)* 1e-3;

weightRange = [0.5 1 1.5 2 2.5 3];


FVFRange = 40;

T2myelRange = 4*1e-3;
T2outRange = 20* 1e-3;

weightRange = 0.5;


nb_TE = 12;
nb_orientations = 6;

noise = 0.01;
nb_replic = 1;

it = 0;
experience_name = 'Porcine2'

options.include_theta = 0;

options.coordinate.classic_polar = 0;
options.coordinate.classic_cartesian = 0;
options.coordinate.polyfit_polar = 0;
options.coordinate.polyfit_cartesian = 1;
options.coordinate.polyfit_cartesian_demean = 0;

for FVF = FVFRange
    for num = 1 : nb_replic
        it = it + 1;
        suffix = ['_train' num2str(num)];
        output_folder = [dico_folder '/dictionary_part/FVF' num2str(FVF) '_N400_train' num2str(num)];
        mkdir(output_folder)
        
        FVF_folder = [signal_folder 'FVF' num2str(FVF) '_N400' suffix '/'];
        signal_path = [FVF_folder 'Signal_FVF' num2str(FVF) suffix '.mat'];
        
        createDictionaryPartSignalMultiModalities(signal_path, output_folder, experience_name, T2myelRange, T2outRange, weightRange, nb_TE, noise, FVF, nb_orientations, num, options);
%         job{it} = qsubfeval(@createDictionaryPartSignalMultiModalities_selective_TE, signal_path, output_folder, T2myelRange, T2outRange, weightRange, nb_TE, noise, FVF, nb_orientations, num, options, 'memreq',  req_mem,  'timreq',  req_etime);
    end
end






