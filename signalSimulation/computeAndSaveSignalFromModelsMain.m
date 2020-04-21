% clear
close all
cd /project/3015069.04/temp/Jobs

model_inputfolder = '/project/3015069.04/WM_Models/N400/';
signal_outputfolder = '/project/3015069.04/signal_components/single_orientation/BrainSample2/';

rotation_folder = '/project/3015069.04/data/rotations/';
req_mem   = 3e9;
req_etime = 2000;

%%%%%%%%%%%% Need to be check
% load([rotation_folder 'theorically_good_16_rotations.mat']);
load([rotation_folder '20_fiber_orientations_rotations.mat']);
% load([rotation_folder 'BrainSample2_rotations_ref_2_orientations.mat'])

% dict_params.TE = linspace(1.8,94.6,30)*1e-3;
% dict_params.TE = linspace(1.7,35.25,12)*1e-3;
% dict_params.TE = linspace(2.4,50.7,12)*1e-3;
% dict_params.TE = (1.8:3.09:54.33)*1e-3;
dict_params.TE = linspace(2.15,35.7,12)*1e-3; % Brain Sample 2 docs

%%%%%%%%%%%%
dict_params.B0 = 3;

dict_params.myelin.xiRange = -0.2:0.1:0.2;
% dict_params.myelin.xaRange = -0.2:0.1:0.2;
dict_params.myelin.xaRange = -0.1;

dict_params.gRatioRange = 50:5:85;

dict_params.fiber_directions = fiber_directions;

dict_params.rotations = 0;
if dict_params.rotations
    dict_params.sphere_rotations = rotations;
end

dict_params.dispersion = 0;
if dict_params.dispersion
    dict_params.kappa_list = [inf 18 9 5.5 3.5];
    dict_params.dispersion_list = [0 0.1 0.2 0.3 0.4];
end

dict_params.no_mask_tensor_map = 0;

it = 0;
nb_replica = 8;

for k = 1:nb_replica
    suffix = ['_train' num2str(k)];
    for FVF = 10 : 10 : 80
        display(['FVF: ' num2str(FVF) ', nb replic: ' num2str(k)]);
        tic()
        filename = [signal_outputfolder 'FVF' num2str(FVF) '_N400' suffix '/' 'Signal_FVF' num2str(FVF) suffix '.mat'];
        out = computeAndSaveSignalFromModels(FVF, suffix, model_inputfolder, signal_outputfolder, dict_params);
        if ~isfile(filename)
            it = it +1
%             job{it} = qsubfeval(@computeAndSaveSignalFromModels, FVF, suffix, model_inputfolder, signal_outputfolder, dict_params, 'memreq',  req_mem,  'timreq',  req_etime);
        end
       
        toc()
    end
end
