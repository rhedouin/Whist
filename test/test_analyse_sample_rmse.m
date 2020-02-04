% Test analyseSampleRmse
cd /project/3015069.04/Jobs

req_mem   = 10e9;
req_memoverhead   = 3e9;
req_etime = 40000;

ses = 'ses-mri01';

parameter_list = {'FVF', 'gRatio', 'xi', 'T2out', 'T2myel', 'weight'};

dict_folder = '/project/3015069.04/dictionaries/multi_orientations/Porcine-1/';
dict_name = 'SignalWithNoise0_1rep_9orientations_Porcine1_fix_xa_polyfit_cartesian_without_theta.h5py';
dict_path = [dict_folder dict_name];

data_folder = ['/project/3015069.01/derived/Porcine-1/' ses '/concatenate_signals/individual_samples/'];
    
nbSample = 1;

it = 0 ;
for kSample = 1:nbSample
    it = it+1
    data_path = [data_folder 'Porcine-1_' ses '_small_sample_' num2str(kSample) '_concatenate_signal_polyfit_cartesian_all_orientations_without_theta_2_ref.nii.gz'];
    mask_path = [data_folder 'Porcine-1_' ses '_small_sample_' num2str(kSample) '_mask.nii.gz'];

    result_path = [data_folder 'Porcine-1_' ses '_small_sample_' num2str(kSample) '_result.txt'];
    options.save_results = result_path;
    
    parameter_map_folder = [data_folder 'parameter_maps/rmse_polyfit_cartesian/'];
    options.save_parameter_map_folder = parameter_map_folder;
    options.save_parameter_map_prefix = ['Porcine-1_' ses '_small_sample_' num2str(kSample)];

%     job{it} = qsubfeval(@analyseSampleRmse, dict_path, data_path, mask_path, parameter_list, options, 'memreq',  req_mem,  'timreq',  req_etime, 'memoverhead',req_memoverhead);
   out = analyseSampleRmse(dict_path, data_path, mask_path, parameter_list, options);

end

% analyseSampleRmse(dict_path, data_path, mask_path, parameter_list, options)