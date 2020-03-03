clear
close all
cd /project/3015069.04/temp/Jobs

model_inputfolder = '/project/3015069.04/WM_Models/N400/';
signal_outputfolder = '/project/3015069.04/signal_components/multi_orientations/Porcine2/fix_xa_large_FVF_single_direction_tensor_mask/';

rotation_folder = '/project/3015069.04/data/rotations/';
req_mem   = 2e9;
req_etime = 2000;

%%%%%%%%%%%% Need to be check
load([rotation_folder 'Porcine-2_rotation_ref_2_orientations.mat']);
load([rotation_folder '20_fiber_orientations_rotations.mat']);

dict_params.TE = linspace(2.44,50.73,12)*1e-3;
%%%%%%%%%%%%

dict_params.B0 = 3;

dict_params.myelin.xiRange = [-0.2, 0.1, 0, 0.1, 0.2];
% dict_params.myelin.xaRange = -0.2:0.1:0.2;
dict_params.myelin.xaRange = -0.1;

dict_params.gRatioRange = 50:5:85;

fiber_directions = [1, 0, 0];

dict_params.fiber_directions = fiber_directions;

dict_params.rotations = 1;
dict_params.sphere_rotations = rotations;

dict_params.dispersion = 0;
dict_params.kappa_list = [inf 18 9 5.5 3.5];
dict_params.dispersion_list = [0 0.1 0.2 0.3 0.4];

it = 0;
nb_replica = 8;

for k = 1:nb_replica
    it = it + 1;
    suffix = ['_train' num2str(k)];
%     for FVF = 40 : 5 : 85
    for FVF = 10 : 10 : 80

%     for FVF = 10 : 10 : 30

        tic()
%         out = computeAndSaveSignalFromModels(FVF, suffix, model_inputfolder, signal_outputfolder, dict_params);
        job{it} = qsubfeval(@computeAndSaveSignalFromModels, FVF, suffix, model_inputfolder, signal_outputfolder, dict_params, 'memreq',  req_mem,  'timreq',  req_etime);
       
        toc()
    end
end
