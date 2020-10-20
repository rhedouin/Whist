% createDictionaryCompletePipeline


cd /project/3015069.04/temp/Jobs
project_folder = '/project/3015069.04/';

model_folder = [project_folder 'WM_Models/N400/'];

experience_name = 'BrainSample2NewRange';

signal_component_folder = [project_folder 'signal_components/multi_orientations/' experience_name '/'];
dict_folder = [project_folder  'dictionaries/multi_orientations/' experience_name '/'];

mkdir(signal_component_folder)

rotation_folder = '/project/3015069.04/data/rotations/';
req_mem1   = 5e9;
req_mem2   = 60e9;

req_etime1 = 5000;
req_etime2 = 5000;

%%%%%%%%%%%% Need to be set one particular set of rotations, the 20 fiber
%%%%%%%%%%%% orientations are evenly spread on the sphere
load([rotation_folder 'BrainSample2_rotations_ref_2_orientations.mat']);
load([rotation_folder '20_fiber_orientations_rotations.mat']);

%%%%%%%%%%%% Signal component parameters
dict_params.TE = linspace(1.7,35.25,12)*1e-3;
nb_TE = length(dict_params.TE);

dict_params.B0 = 3;

dict_params.myelin.xiRange = -0.2:0.1:0.2;
dict_params.myelin.xaRange = -0.1;

dict_params.gRatioRange = 50:5:85;

dict_params.fiber_directions = fiber_directions;

dict_params.rotations = 1;
if dict_params.rotations == 0
    nb_rotations = 1;
else
   dict_params.sphere_rotations = rotations;
   nb_rotations = size(rotations, 3);
end


dict_params.dispersion = 0;
dict_params.kappa_list = [inf 18 9 5.5 3.5];
dict_params.dispersion_list = [0 0.1 0.2 0.3 0.4];


%%%%%%%%%%%% Final signal parameters
FVFRange = (10 : 10 : 80);

T2MyelinRange = [4 8 12 16 20]*1e-3;
T2IntraExtraAxonalRange = (15 : 10 : 95)* 1e-3;
T2ExtraAxonalRange = (15 : 10 : 95)* 1e-3;

weightRange = [0.5 1 1.5 2 2.5 3];

noise = 0.04;

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
                
        dict_params.FVF_round = FVF;

        tic()
        model_FVF_folder = [model_folder 'FVF' num2str(FVF) '_N400' model_suffix '/'];

        signal_component_FVF_folder = [signal_component_folder 'FVF' num2str(FVF) '_N400' model_suffix '/'];
        mkdir(signal_component_FVF_folder)

        signal_component_FVF_path = [signal_component_FVF_folder 'Signal_FVF' num2str(FVF) model_suffix '.mat'];
        
        dict_FVF_folder = [dict_folder 'dictionary_part/FVF' num2str(FVF) '_N400_train' num2str(k) '/'];
        mkdir(dict_FVF_folder)

%          computeCompartmentSignalFromModels(model_folder, signal_component_FVF_path, model_suffix, dict_params);
        jobId1 = qsubfeval(@computeCompartmentSignalFromModels, model_folder, signal_component_FVF_path, model_suffix, dict_params, 'memreq',  req_mem1,  'timreq',  req_etime1);
        jobIdAll{end+1} = jobId1;

%         computeTotalSignalDictionaryPart(signal_component_FVF_path, dict_FVF_folder, experience_name, T2MyelinRange, T2IntraExtraAxonalRange, weightRange, nb_TE, noise, FVF, nb_rotations, k, options);

        jobId2 = qsubfeval(@computeTotalSignalDictionaryPart, signal_component_FVF_path, dict_FVF_folder, experience_name, T2MyelinRange, T2IntraExtraAxonalRange, weightRange, nb_TE, noise, FVF, nb_rotations, k, options, 'memreq',  req_mem1,  'timreq',  req_etime1, 'waitfor', jobId1);
        jobIdAll{end+1} = jobId2;  
  
        toc()
    end
end

prefix_name = ['SignalWithNoise' num2str(100*noise)];

qsubfeval(@concatenateDictionaryParts, dict_folder, FVFRange, prefix_name, nb_replica, nb_rotations, nb_TE, experience_name, dict_suffix, 'memreq',  req_mem2,  'timreq',  req_etime2, 'waitfor', jobIdAll)



