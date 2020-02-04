close all
clear

dict_folder = '/project/3015069.04/dictionaries/multi_orientations/BrainSample2/';
data_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals/';
data_path = [data_folder 'BrainSample-2_ses-03_all_orientations_fa-20_concatenate_signal_2_orientation-4.nii.gz'];

% load dictionary
tic()
dict_name = 'SignalWithNoise0_1rep_with_theta_9orientations_BrainSample2_fix_xa.h5py';
dict_path = [dict_folder dict_name];
dict_signal_values = h5read(dict_path, '/SignalValues');
dict_signal_shape = size(dict_signal_values);
signal_length = dict_signal_shape(1);
dict_shape = dict_signal_shape(2:end);
dict_length = prod(dict_shape);

FVF_grid =    h5read(dict_path, '/FVFValues');
gRatio_grid = h5read(dict_path, '/gRatioValues');
xi_grid =     h5read(dict_path, '/xiValues');
T2out_grid =    h5read(dict_path, '/T2outValues');
T2myel_grid =    h5read(dict_path, '/T2myelValues');
weight_grid =    h5read(dict_path, '/weightValues');
direction_grid = h5read(dict_path, '/directionsValues');

dict_signal_reshape = reshape(dict_signal_values, signal_length, dict_length);
FVF_reshape = FVF_grid(:);
gRatio_reshape = gRatio_grid(:);
xi_reshape = xi_grid(:);
T2out_reshape = T2out_grid(:);
T2myel_reshape = T2myel_grid(:);
weight_reshape = weight_grid(:);


% load log scale dictionary

tic()
dict_log_name = 'SignalWithNoise0_1rep_with_theta_9orientations_BrainSample2_fix_xa_log_scale.h5py';
dict_log_path = [dict_folder dict_log_name];
dict_log_signal_values = h5read(dict_log_path, '/SignalValues');
dict_log_signal_shape = size(dict_log_signal_values);

signal_log_length = dict_log_signal_shape(1);
dict_log_shape = dict_log_signal_shape(2:end);
dict_log_length = prod(dict_log_shape);

FVF_log_grid =    h5read(dict_log_path, '/FVFValues');
gRatio_log_grid = h5read(dict_log_path, '/gRatioValues');
xi_log_grid =     h5read(dict_log_path, '/xiValues');
T2out_log_grid =    h5read(dict_log_path, '/T2outValues');
T2myel_log_grid =    h5read(dict_log_path, '/T2myelValues');
weight_log_grid =    h5read(dict_log_path, '/weightValues');
direction_log_grid = h5read(dict_log_path, '/directionsValues');

dict_log_signal_reshape = reshape(dict_log_signal_values, signal_log_length, dict_log_length);
FVF_log_reshape = FVF_log_grid(:);
gRatio_log_reshape = gRatio_log_grid(:);
xi_log_reshape = xi_log_grid(:);
T2out_log_reshape = T2out_log_grid(:);
T2myel_log_reshape = T2myel_log_grid(:);
weight_log_reshape = weight_log_grid(:);

% Load data 
data_signal_values = load_nii_img_only(data_path);
temp_shape = size(data_signal_values);
data_shape = temp_shape(1:3);
data_signal_reshape = permute(data_signal_values,[4 1 2 3]);

% data_length = prod(data_shape);

% mask_path = [data_folder '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_gre_orientation-4_brain_mask.nii.gz'];
% mask = load_nii_img_only(mask_path);

% direction_reshape = reshape(direction_grid, [3 dict_length]);
toc()
% data_reshape = reshape(data_reshape, signal_length, data_length);
nb_map = 1;

CSF_list_pixel = {[53, 69, 65], [74, 70, 65], [54, 69, 61], [72, 67, 61]}

for k = 1:length(CSF_list_pixel)
    data_signal_selection = data_signal_reshape(:, CSF_list_pixel{k}(1), CSF_list_pixel{k}(2), CSF_list_pixel{k}(3));
    
    data_signal_replicate = repmat(data_signal_selection,1,dict_length);
    mse = sum((data_signal_replicate - dict_signal_reshape).^2, 1);
    
    [mse_max, mse_sort] = sort(mse);
    
    figure
    plot(data_signal_selection, 'k')
    hold on
    
    for nb = 1:nb_map
        best_dict_match = dict_signal_reshape(:, mse_sort(nb))
        
        FVF_CSF(k,nb) = FVF_reshape(mse_sort(nb));
        gRatio_CSF(k,nb) = gRatio_reshape(mse_sort(nb));
        xi_CSF(k,nb) = xi_reshape(mse_sort(nb));
        T2out_CSF(k,nb) = T2out_reshape(mse_sort(nb));
        T2myel_CSF(k,nb) = T2myel_reshape(mse_sort(nb));
        weight_CSF(k,nb) = weight_reshape(mse_sort(nb));
        
        plot(best_dict_match, 'b')        
        
    end
    
    data_log_signal_replicate = repmat(data_signal_selection,1,dict_log_length);
    mse_log = sum((data_log_signal_replicate - dict_log_signal_reshape).^2, 1);
    
    [mse_log_max, mse_log_sort] = sort(mse_log);
       
    for nb = 1:nb_map
        best_dict_log_match = dict_log_signal_reshape(:, mse_log_sort(nb))
        
        FVF_log_CSF(k,nb) = FVF_log_reshape(mse_log_sort(nb));
        gRatio_log_CSF(k,nb) = gRatio_log_reshape(mse_log_sort(nb));
        xi_log_CSF(k,nb) = xi_log_reshape(mse_log_sort(nb));
        T2out_log_CSF(k,nb) = T2out_log_reshape(mse_log_sort(nb));
        T2myel_log_CSF(k,nb) = T2myel_log_reshape(mse_log_sort(nb));
        weight_log_CSF(k,nb) = weight_log_reshape(mse_log_sort(nb));
        
        plot(best_dict_log_match, 'r')
        title({['CSF classic, FVF : ' num2str(FVF_CSF(k, nb)) ' g ratio : ' num2str(gRatio_CSF(k, nb)) ...
            ' T2out: ' num2str(T2out_CSF(k, nb)) ' T2myel : ' num2str(T2myel_CSF(k, nb)) ' xi : ' num2str(xi_CSF(k, nb)) ...
            ' weight : ' num2str(weight_CSF(k, nb)) ' mse sort : ' num2str(mse_max(nb))], ...
            ['log, FVF : ' num2str(FVF_log_CSF(k, nb)) ' g ratio : ' num2str(gRatio_log_CSF(k, nb)) ...
            ' T2out: ' num2str(T2out_log_CSF(k, nb)) ' T2myel : ' num2str(T2myel_log_CSF(k, nb)) ' xi : ' num2str(xi_log_CSF(k, nb)) ...
            ' weight : ' num2str(weight_log_CSF(k, nb)) ' mse sort : ' num2str(mse_log_max(nb))]})           
        
        legend('data', 'classic scale', 'log scale')
    end
end

CC_list_pixel = {[74, 75, 61], [54, 57, 61], [75, 75, 65], [53, 77, 65]}

for k = 1:length(CC_list_pixel)
    data_signal_selection = data_signal_reshape(:, CC_list_pixel{k}(1), CC_list_pixel{k}(2), CC_list_pixel{k}(3));
    
    data_signal_replicate = repmat(data_signal_selection,1,dict_length);
    mse = sum((data_signal_replicate - dict_signal_reshape).^2, 1);
    
    [mse_max, mse_sort] = sort(mse);
    
    nb = 1;     
    
    figure
    plot(data_signal_selection, 'k')
    hold on
        
    for nb = 1:nb_map
        best_dict_match = dict_signal_reshape(:, mse_sort(nb))
        
        FVF_CC(k,nb) = FVF_reshape(mse_sort(nb));
        gRatio_CC(k,nb) = gRatio_reshape(mse_sort(nb));
        xi_CC(k,nb) = xi_reshape(mse_sort(nb));
        T2out_CC(k,nb) = T2out_reshape(mse_sort(nb));
        T2myel_CC(k,nb) = T2myel_reshape(mse_sort(nb));
        weight_CC(k,nb) = weight_reshape(mse_sort(nb));
        
        plot(best_dict_match, 'b')
        
    end
    
    data_log_signal_replicate = repmat(data_signal_selection,1,dict_log_length);
    mse_log = sum((data_log_signal_replicate - dict_log_signal_reshape).^2, 1);
    
    [mse_log_max, mse_log_sort] = sort(mse_log);
       
    for nb = 1:nb_map
        best_dict_log_match = dict_log_signal_reshape(:, mse_log_sort(nb))
        
        FVF_log_CC(k,nb) = FVF_log_reshape(mse_log_sort(nb));
        gRatio_log_CC(k,nb) = gRatio_log_reshape(mse_log_sort(nb));
        xi_log_CC(k,nb) = xi_log_reshape(mse_log_sort(nb));
        T2out_log_CC(k,nb) = T2out_log_reshape(mse_log_sort(nb));
        T2myel_log_CC(k,nb) = T2myel_log_reshape(mse_log_sort(nb));
        weight_log_CC(k,nb) = weight_log_reshape(mse_log_sort(nb));
        
        plot(best_dict_log_match, 'r')
        title({['CC classic, FVF : ' num2str(FVF_CC(k, nb)) ' g ratio : ' num2str(gRatio_CC(k, nb)) ...
            ' T2out: ' num2str(T2out_CC(k, nb)) ' T2myel : ' num2str(T2myel_CC(k, nb)) ' xi : ' num2str(xi_CC(k, nb)) ...
            ' weight : ' num2str(weight_CC(k, nb)) ' mse sort : ' num2str(mse_max(nb))], ...
            ['log, FVF : ' num2str(FVF_log_CC(k, nb)) ' g ratio : ' num2str(gRatio_log_CC(k, nb)) ...
            ' T2out: ' num2str(T2out_log_CC(k, nb)) ' T2myel : ' num2str(T2myel_log_CC(k, nb)) ' xi : ' num2str(xi_log_CC(k, nb)) ...
            ' weight : ' num2str(weight_log_CC(k, nb)) ' mse sort : ' num2str(mse_log_max(nb))]})   
        
        legend('data', 'classic scale', 'log scale') 
    end
end

