% createDictionaryCompletePipeline


cd /project/3015069.04/temp/Jobs
project_folder = '/project/3015069.04/';

model_folder = [project_folder 'WM_Models/N400/'];
signal_component_folder = [project_folder 'signal_components/multi_orientations/BS-3/'];
dict_folder = [project_folder  'dictionaries/multi_orientations/BS-3/'];

mkdir(signal_component_folder)

rotation_folder = '/project/3015069.04/data/rotations/';
req_mem1   = 6e9;
req_mem2   = 50e9;

req_etime = 10000;

%%%%%%%%%%%% Need to be check
load([rotation_folder 'BS3_highres_rotations_ref_2_orientations.mat']);
load([rotation_folder '20_fiber_orientations_rotations.mat']);

%%%%%%%%%%%% Signal component parameters
dict_params.TE = linspace(2.12,54.31,18)*1e-3;

dict_params.B0 = 3;

dict_params.myelin.xiRange = -0.2:0.1:0.2;
dict_params.myelin.xaRange = -0.1;

dict_params.gRatioRange = 50:5:85;

dict_params.fiber_directions = fiber_directions;

dict_params.rotations = 1;
dict_params.sphere_rotations = rotations;

dict_params.dispersion = 0;
dict_params.kappa_list = [inf 18 9 5.5 3.5];
dict_params.dispersion_list = [0 0.1 0.2 0.3 0.4];


%%%%%%%%%%%% Final signal parameters
FVFRange = (10 : 10 : 80);

T2MyelinRange = [4 8 12 16 20]*1e-3;
T2IntraExtraAxonalRange = (20 : 20 : 100)* 1e-3;
T2ExtraAxonalRange = (20 : 20 : 100)* 1e-3;

weightRange = [0.5 1 1.5 2 2.5 3];

nb_TE = 18;
nb_rotations = 10;

noise = 0.01;
nb_replic = 8;

experience_name = 'BS-3';

options.include_theta = 1;

options.coordinate.polyfit_cartesian = 1;
options.coordinate.polyfit_cartesian_demean = 0;
options.coordinate.polyfit_polar = 0;
options.coordinate.classic_cartesian = 0;
options.coordinate.classic_polar = 0;

dict_suffix = '_polyfit_cartesian_with_theta';

it = 0;
nb_replica = 8;
jobIdAll = {};

for k = 1:nb_replica
    it = it + 1;
    model_suffix = ['_train' num2str(k)];
    for FVF = FVFRange
                
        tic()
        model_FVF_folder = [model_folder 'FVF' num2str(FVF) '_N400' model_suffix '/'];

        signal_component_FVF_folder = [signal_component_folder 'FVF' num2str(FVF) '_N400' model_suffix '/'];
        mkdir(signal_component_FVF_folder)

        signal_component_FVF_path = [signal_component_FVF_folder 'Signal_FVF' num2str(FVF) model_suffix '.mat'];
        
        dict_FVF_folder = [dict_folder 'dictionary_part/FVF' num2str(FVF) '_N400_train' num2str(k) '/'];
        mkdir(dict_FVF_folder)
        
%         computeAndSaveSignalFromModels(FVF, model_suffix, model_FVF_folder, signal_component_FVF_path, dict_params);
%         jobId = qsubfeval(@computeAndSaveSignalFromModels, FVF, model_suffix, model_FVF_folder, signal_component_FVF_path, dict_params, 'memreq',  req_mem,  'timreq',  req_etime);
%         jobIdAll{end+1} = jobId;

%         createDictionaryPartSignalMultiModalities(signal_component_FVF_path, dict_FVF_folder, experience_name, T2MyelinRange, T2IntraExtraAxonalRange, weightRange, nb_TE, noise, FVF, nb_rotations, k, options);
%         jobId2 = qsubfeval(@createDictionaryPartSignalMultiModalities, signal_component_FVF_path, dict_FVF_folder, experience_name, T2MyelinRange, T2IntraExtraAxonalRange, weightRange, nb_TE, noise, FVF, nb_rotations, k, options, 'memreq',  req_mem,  'timreq',  req_etime, 'waitfor', jobId);
%         jobId2 = qsubfeval(@createDictionaryPartSignalMultiModalities, signal_component_FVF_path, dict_FVF_folder, experience_name, T2MyelinRange, T2IntraExtraAxonalRange, weightRange, nb_TE, noise, FVF, nb_rotations, k, options, 'memreq',  req_mem1,  'timreq',  req_etime);
% 
%         jobIdAll{end+1} = jobId2;    
   
        toc()
    end
end

prefix_name = ['SignalWithNoise' num2str(100*noise)];

concatenateDictionaryPartsFunction(dict_folder, FVFRange, prefix_name, nb_replica, nb_rotations, nb_TE, experience_name, dict_suffix)
% qsubfeval(@concatenateDictionaryPartsFunction, dict_folder, FVFRange, prefix_name, nb_replica, nb_rotations, nb_TE, experience_name, dict_suffix, 'memreq',  req_mem2,  'timreq',  req_etime, 'waitfor', jobIdAll)




