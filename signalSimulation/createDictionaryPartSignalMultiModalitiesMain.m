
clear
% close all

req_mem   = 5e9;
req_etime = 6000;

base_folder = '/project/3015069.04/';
job_folder = [base_folder 'temp/Jobs'];
signal_folder = [base_folder 'signal_components/multi_orientations/BrainSample2/'];
dico_folder = [base_folder 'dictionaries/multi_orientations/BrainSample2/'];

cd(job_folder)

FVFRange = (10 : 10 : 80);

T2MyelinRange = [4 8 12 16 20]*1e-3;
T2IntraExtraAxonalRange = (20 : 20 : 100)* 1e-3;
T2ExtraAxonalRange = (20 : 20 : 100)* 1e-3;

weightRange = [0.5 1 1.5 2 2.5 3];

nb_TE = 12;
nb_orientations = 9;

noise = 0;
nb_replic = 8;

it = 0;
experience_name = 'BrainSample2';

options.include_theta = 1;

options.coordinate.classic_polar = 0;
options.coordinate.classic_cartesian = 0;
options.coordinate.polyfit_polar = 0;
options.coordinate.polyfit_cartesian = 1;
options.coordinate.polyfit_cartesian_demean = 0;

noise_list = [0.01];

replic_list = [1, 2, 3, 4, 5, 6, 7, 8];

coordinate = 'polyfit_cartesian';

for noise = noise_list
    for FVF = FVFRange
        for k = 1 : length(replic_list)
            num = replic_list(k);
            it = it + 1
            suffix = ['_train' num2str(num)];
            output_folder = [dico_folder 'dictionary_part/FVF' num2str(FVF) '_N400_train' num2str(num) '/'];
            mkdir(output_folder)
            
            FVF_folder = [signal_folder 'FVF' num2str(FVF) '_N400' suffix '/'];
            signal_path = [FVF_folder 'Signal_FVF' num2str(FVF) suffix '.mat'];
            %                 createDictionaryPartSignalMultiModalities(signal_path, output_folder, experience_name, T2MyelinRange, T2IntraExtraAxonalRange, weightRange, nb_TE, noise, FVF, nb_orientations, num, options);
            
            if noise == 0.005
                prefix_name = 'SignalWithNoise05';
            else
                prefix_name = ['SignalWithNoise' num2str(100*noise)];
            end
            
            suffix_theta = 'with_theta';
            
            base_name = [prefix_name '_FVF' num2str(FVF) '_replic' num2str(num) '_' num2str(nb_orientations) ...
                'orientations_' num2str(nb_TE)  'TE_' experience_name '_fix_xa_' coordinate '_'  suffix_theta];
            signal_name = [base_name '.h5py'];
            
            filename =[output_folder  signal_name];
            
            if ~isfile(filename)
                filename
                tic()
%                   createDictionaryPartSignalMultiModalities(signal_path, output_folder, experience_name, T2MyelinRange, T2IntraExtraAxonalRange, weightRange, nb_TE, noise, FVF, nb_orientations, num, options);
                job{it} = qsubfeval(@createDictionaryPartSignalMultiModalities, signal_path, output_folder, experience_name, T2MyelinRange, T2IntraExtraAxonalRange, weightRange, nb_TE, noise, FVF, nb_orientations, num, options, 'memreq',  req_mem,  'timreq',  req_etime);
                toc()
            end
        end
    end
end




