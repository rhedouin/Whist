% close all
clear

dict_folder = '/project/3015069.04/dictionaries/multi_orientations/Porcine2/fix_xa_large_FVF_20_directions/';
data_folder = '/project/3015069.01/derived/Porcine-2/ses-mri01/concatenate_signals/individual_samples/';
rmse_folder = [data_folder 'parameter_maps/rmse/'];

parameter_list = {'FVF', 'gRatio', 'xiMyelin', 'T2IntraExtraAxonal', 'T2Myelin', 'weight'};
lPara = length(parameter_list);
%%%%%% Load dictionary signal values
tic()
dict_name = 'SignalWithNoise0_1rep_6orientations_12TE_Porcine2_fix_xa_large_FVF_20_directions_polyfit_cartesian_with_theta.h5py';
dict_path = [dict_folder dict_name];
dict_signal_values = h5read(dict_path, '/SignalValues');
dict_signal_shape = size(dict_signal_values);
signal_length = dict_signal_shape(1);
dict_shape = dict_signal_shape(2:end);
dict_length = prod(dict_shape);
dict_signal_reshape = reshape(dict_signal_values, signal_length, dict_length);

normalize = 1;
if normalize == 1
    max_signal = max(abs(dict_signal_reshape), [], 2);
else
    max_signal = ones(size(dict_signal_reshape, 1), 1);
end

dict_signal_reshape_norm = dict_signal_reshape./max_signal;

%%%%%% Load parameter grid
for kPara = 1:lPara
    parameter = parameter_list{kPara};
    parameter_grid.(parameter) = h5read(dict_path, ['/' parameter 'Values']);
    parameter_reshape.(parameter) = parameter_grid.(parameter)(:);
end

%%%%%%%%%%%% Load data
for kSample = 1:7
    
    nSample = kSample;    
    experience_name = 'Porcine-2_ses-mri01';
    
    % with theta
    display('With theta')
    sample_name = [experience_name '_small_sample_' num2str(nSample) '_concatenate_signal_all_orientations_polyfit_cartesian_with_theta.nii.gz'];
    data_path = [data_folder sample_name];
    
    mask_name = [experience_name '_small_sample_' num2str(nSample) '_mask.nii.gz'];
    mask_path = [data_folder mask_name];
    se = strel('cube',1)
    
    mask = single(load_nii_img_only(mask_path));
    % mask_erode = imerode(mask, se);
    mask_erode = mask;
%     mask_erode = zeros(size(mask));
    
    % Load data
    data_signal_nii = load_untouch_nii(data_path);
    data_signal_values = data_signal_nii.img;
    
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
                    
%                     subplot(221)
%                     plot(data_signal_selection)
                    
                    data_avg = data_avg + data_signal_selection;
                    hold on
                    
                    % Compute rmse for each pixel
                    data_signal_selection_norm = data_signal_selection./max_signal;
                    data_signal_replicate_norm = repmat(data_signal_selection_norm,1,dict_length);
                    mse = sum((data_signal_replicate_norm - dict_signal_reshape_norm).^2, 1);
                    
                    [mse_max, mse_sort] = sort(mse);
                    
                    for kPara = 1:lPara
                        parameter = parameter_list{kPara};
                        parameter_map.(parameter)(k, l, m) = parameter_reshape.(parameter)(mse_sort(1));
                        parameter_map_avg.(parameter) = parameter_map_avg.(parameter) + parameter_map.(parameter)(k, l, m);
                    end
                    
                    it = it + 1
                    
                end
            end
        end
    end
    toc()
    for kPara = 1:lPara
        parameter = parameter_list{kPara};
        parameter_map_avg.(parameter) = parameter_map_avg.(parameter) / it;
        
        title_parameter_map = [];
    end
   
    % Save parameter maps
    for kPara = 1:lPara
        parameter = parameter_list{kPara};
        parameter_map_name = [experience_name '_small_sample_' num2str(nSample) '_' parameter '_rmse.nii.gz'];
        data_signal_nii.img = single(parameter_map.(parameter));
        
        save_untouch_nii(data_signal_nii, [rmse_folder parameter_map_name]);
    end
    
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
%         
%     subplot(212)
%     plot(data_avg)
%     hold on
%     
%     plot(best_dict_match)
%     legend('data', 'dict best match')
%     
%     title(title_data_avg)
    
end

