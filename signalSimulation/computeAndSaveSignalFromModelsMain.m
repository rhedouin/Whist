clear
close all
cd /project/3015069.04/temp/Jobs

model_inputfolder = '/project/3015069.04/WM_Models/N400/';
signal_outputfolder = '/project/3015069.04/signal_components/multi_orientations/Porcine2/';

rotation_folder = '/project/3015069.04/data/rotations/';
req_mem   = 3e9;
req_etime = 1000;

%%%%%%%%%%%% Need to be check
load([rotation_folder 'Porcine-2_rotation_ref_2_orientations.mat']);
% load([rotation_folder '20_fiber_orientations_rotations.mat']);

options.time = linspace(2.44,50.73,12)*1e-3;
%%%%%%%%%%%%

options.B0 = 3;

options.xiRange = -0.2:0.1:0.2;
options.xaRange = -0.1;
options.gRatioRange = 50:5:85;

fiber_directions = [1, 0, 0];

options.fiber_directions = fiber_directions;

options.rotations = 1;
options.sphere_rotations = rotations;

options.dispersion = 0;
options.kappa_list = [inf 18 9 5.5 3.5];
options.dispersion_list = [0 0.1 0.2 0.3 0.4];

it = 0;
for k = 1:8
    it = it + 1;
    suffix = ['_train' num2str(k)];
    for FVF = 40 : 5 : 85
%     for FVF = 10 : 10 : 30

        tic()
        job{it} = qsubfeval(@computeAndSaveSignalFromModels, FVF, suffix, model_inputfolder, signal_outputfolder, options, 'memreq',  req_mem,  'timreq',  req_etime);
%         out = computeAndSaveSignalFromModels(FVF, suffix, model_inputfolder, signal_outputfolder, options);
        
        toc()
    end
end
