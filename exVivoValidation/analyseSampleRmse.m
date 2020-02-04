function out = analyseSampleRmse(dict_path, data_path, mask_path, parameter_list, options)
    
lPara = length(parameter_list);
%%%%%% Load dictionary signal values
tic()
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

%%%%%% Load data and mask
data_signal_values_nii = load_untouch_nii(data_path);
data_signal_values = data_signal_values_nii.img;

mask = single(load_nii_img_only(mask_path));
nb_pixel = sum(mask, 'all');

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
display(['nb of pixel to process : ' num2str(nb_pixel)]);
for k=1:dims(1)
    for l=1:dims(2)
        for m=1:dims(3)
            if mask(k, l, m) == 1
                
                data_signal_selection = data_signal_reshape(:, k, l, m);                
                data_avg = data_avg + data_signal_selection;

                % Compute rmse for each pixel
                data_signal_replicate = repmat(data_signal_selection,1,dict_length);
                mse = sum((data_signal_replicate - dict_signal_reshape).^2, 1);
                
                [~, mse_sort] = sort(mse);
                
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
    parameter_map_avg.(parameter) = parameter_map_avg.(parameter) / nb_pixel;
end

% Compute rmse for average signal
data_avg = data_avg / nb_pixel;
data_signal_selection = data_avg;

data_signal_replicate = repmat(data_signal_selection,1,dict_length);
mse = sum((data_signal_replicate - dict_signal_reshape).^2, 1);

[~, mse_sort] = sort(mse);

if isfield(options, 'save_parameter_map_folder')
    mkdir(options.save_parameter_map_folder);
    for kPara = 1:lPara
        parameter = parameter_list{kPara};
        output_path = [options.save_parameter_map_folder options.save_parameter_map_prefix '_' parameter '_rmse_polyfit.nii.gz'];
        
        data_signal_values_nii.img = parameter_map.(parameter);
        save_untouch_nii(data_signal_values_nii, output_path);
    end
end

if isfield(options, 'save_results')     
    fileID = fopen(options.save_results, 'w');
    
    fprintf(fileID, '%s %s %s \n', '    ', 'Parameter avg', 'Data avg'); 
    for kPara = 1:lPara
        parameter = parameter_list{kPara};
        parameter_data_avg.(parameter) = parameter_reshape.(parameter)(mse_sort(1));
        fprintf(fileID, '%s %f %f \n', parameter, parameter_map_avg.(parameter), parameter_data_avg.(parameter));
    end
    fclose(fileID);
  
end
out = 0;
end
