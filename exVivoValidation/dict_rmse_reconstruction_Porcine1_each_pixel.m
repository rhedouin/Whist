close all
clear

dict_folder = '/project/3015069.04/dictionaries/multi_orientations/Porcine-1/';
data_folder = '/project/3015069.01/derived/Porcine-1/ses-mri01/concatenate_signals/individual_samples/';
    
parameter_list = {'FVF', 'gRatio', 'xi', 'T2out', 'T2myel', 'weight'};
lPara = length(parameter_list);
%%%%%% Load dictionary signal values
tic()
dict_name = 'SignalWithNoise0_1rep_9orientations_Porcine1_fix_xa_polyfit_cartesian_without_theta.h5py';
dict_path = [dict_folder dict_name];
dict_signal_values = h5read(dict_path, '/SignalValues');
dict_signal_shape = size(dict_signal_values);
signal_length = dict_signal_shape(1);
dict_shape = dict_signal_shape(2:end);
dict_length = prod(dict_shape);
dict_signal_reshape = reshape(dict_signal_values, signal_length, dict_length);

%%%%%% Load parameter grid
for kPara = 1:lPara
    parameter = parameter_list{kPara};
    parameter_grid.(parameter) = h5read(dict_path, ['/' parameter 'Values']);
    parameter_reshape.(parameter) = parameter_grid.(parameter)(:);
end

%%%%%%%%%%%% Load data
for kSample = 5:5
figure(kSample)

nSample = kSample;

% without theta
display('Without theta')
data_path = [data_folder 'Porcine-1_ses-mri01_small_sample_' num2str(nSample) '_concatenate_signal_polyfit_cartesian_all_orientations_without_theta_2_ref.nii.gz'];
mask_path = [data_folder 'Porcine-1_ses-mri01_small_sample_' num2str(nSample) '_mask.nii.gz'];

se = strel('cube',1)

mask = single(load_nii_img_only(mask_path));
% mask_erode = imerode(mask, se);
mask_erode = mask;

% Load data 
data_signal_values = load_nii_img_only(data_path);

temp_shape = size(data_signal_values);
dims = temp_shape(1:3);
data_signal_reshape = permute(data_signal_values,[4 1 2 3]);


for kPara = 1:lPara
    parameter = parameter_list{kPara};
    parameter_map.(parameter) = zeros(dims);
    parameter_map_avg.(parameter) = 0;
end

data_avg = 0;
it = 0;
tic()
for k=1:dims(1)
    for l=1:dims(2)
        for m=1:dims(3)
            if mask_erode(k, l, m) == 1
                
                data_signal_selection = data_signal_reshape(:, k, l, m);
                
                subplot(221)
                plot(data_signal_selection)
                
                data_avg = data_avg + data_signal_selection;
                hold on
                
                % Compute rmse for each pixel
%                 data_signal_replicate = repmat(data_signal_selection,1,dict_length);
%                 mse = sum((data_signal_replicate - dict_signal_reshape).^2, 1);
%                 
%                 [mse_max, mse_sort] = sort(mse);
%                 
%                 for kPara = 1:lPara
%                     parameter = parameter_list{kPara};
%                     parameter_map.(parameter)(k, l, m) = parameter_reshape.(parameter)(mse_sort(1));
%                     parameter_map_avg.(parameter) = parameter_map_avg.(parameter) + parameter_map.(parameter)(k, l, m);
%                 end
                
                it = it + 1
                
            end
        end
    end
end
toc()
% 
% for kPara = 1:lPara
%     parameter = parameter_list{kPara};
%     parameter_map_avg.(parameter) = parameter_map_avg.(parameter) / it;
% end

% Compute rmse for average signal
data_avg = data_avg / it;
data_signal_selection = data_avg;

data_signal_replicate = repmat(data_signal_selection,1,dict_length);
mse = sum((data_signal_replicate - dict_signal_reshape).^2, 1);

[mse_max, mse_sort] = sort(mse);

best_dict_match = dict_signal_reshape(:, mse_sort(1));

% title_parameter_avg =  'Parameter average, ';
title_data_avg =  'Data average, ';

for kPara = 1:lPara
    parameter = parameter_list{kPara};
    parameter_data_avg.(parameter) = parameter_reshape.(parameter)(mse_sort(1));
    
%     title_parameter_avg = [title_parameter_avg ', ' parameter ' : ' num2str(parameter_map_avg.(parameter))];
    title_data_avg = [title_data_avg ', ' parameter ' : ' num2str(parameter_data_avg.(parameter))];
end

keyboard;

subplot(212)
plot(data_avg)
hold on

plot(best_dict_match)
legend('data', 'dict best match')

title(title_data_avg)  

subplot(211)
title(title_parameter_avg)  

end

