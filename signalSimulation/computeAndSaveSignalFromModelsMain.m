clear
close all
cd /project/3015069.04/temp/Jobs

model_inputfolder = '/project/3015069.04/WM_Models/N400/';
signal_outputfolder = '/project/3015069.04/signal_components/multi_orientations/theorically_good_15_rotations/';

rotation_folder = '/project/3015069.04/data/rotations/';
req_mem   = 3e9;
req_etime = 10000;

load([rotation_folder 'theorically_good_15_rotations.mat']);
load([rotation_folder '20_fiber_orientations_rotations.mat']);

options.B0 = 3;
options.time = linspace(1.7,35.25,12)*1e-3;
% options.time = (1.7 + 3.05*(0:19))*1e-3;

options.xiRange = -0.2:0.1:0.2;
options.xaRange = -0.1;
options.gRatioRange = 50:5:85;

options.fiber_directions = fiber_directions;

options.rotations = 1;
options.sphere_rotations = rotations;

options.dispersion = 0;
options.kappa_list = [inf 18 9 5.5 3.5];
options.dispersion_list = [0 0.1 0.2 0.3 0.4];

it = 0;
for k = 1:5
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
