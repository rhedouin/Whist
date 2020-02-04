% close all
clear

dict_folder = '/project/3015069.04/dictionaries/multi_orientations/Porcine-1/';
data_folder = '/project/3015069.01/derived/Porcine-1/ses-mri01/concatenate_signals/individual_samples/';
    
for kSample = 5:5
%     figure(kSample)
%     figure

nSample = kSample;

%%%%%%%%%%%%%%%% Without theta
display('Without theta')
data_path = [data_folder 'Porcine-1_ses-mri01_small_sample_' num2str(nSample) '_concatenate_signal_polyfit_cartesian_all_orientations_without_theta_2_ref.nii.gz'];
mask_path = [data_folder 'Porcine-1_ses-mri01_small_sample_' num2str(nSample) '_mask_high_fa.nii.gz'];

se = strel('cube',1)

mask = single(load_nii_img_only(mask_path));
% mask_erode = imerode(mask, se);
mask_erode = mask;

% Load data 
data_signal_values = load_nii_img_only(data_path);

temp_shape = size(data_signal_values);
dims = temp_shape(1:3);
data_signal_reshape = permute(data_signal_values,[4 1 2 3]);
data_avg = 0;
for k=1:dims(1)
    for l=1:dims(2)
        for m=1:dims(3)
            if mask_erode(k, l, m) == 1
                subplot(221)
                plot(data_signal_reshape(:,k,l,m))
                
                data_avg = data_avg + data_signal_reshape(:,k,l,m);
                hold on
            end
        end
    end
end

data_avg = data_avg / sum(mask_erode, 'all');
subplot(222)
plot(data_avg)
hold on

% load dictionary 
tic()
dict_name = 'SignalWithNoise0_1rep_9orientations_Porcine1_fix_xa_polyfit_cartesian_without_theta.h5py';
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

data_signal_selection = data_avg;

data_signal_replicate = repmat(data_signal_selection,1,dict_length);
mse = sum((data_signal_replicate - dict_signal_reshape).^2, 1);

[mse_max, mse_sort] = sort(mse);

nb_map = 1;
for nb = 1:nb_map
    best_dict_match = dict_signal_reshape(:, mse_sort(nb));
    
    FVF(nb) = FVF_reshape(mse_sort(nb));
    gRatio(nb) = gRatio_reshape(mse_sort(nb));
    xi(nb) = xi_reshape(mse_sort(nb));
    T2out(nb) = T2out_reshape(mse_sort(nb));
    T2myel(nb) = T2myel_reshape(mse_sort(nb));
    weight(nb) = weight_reshape(mse_sort(nb));
    
    plot(best_dict_match)
    legend('data', 'dict best match')
    
    title(['FVF : ' num2str(FVF(nb)) ', gRatio : ' num2str(gRatio(nb)) ', xi : ' num2str(xi(nb)) ...
         ', T2out : ' num2str(T2out(nb)) ', T2myel : ' num2str(T2myel(nb)) ', weight : ' num2str(weight(nb))])
    
end


%%%%%%%%%%%%%%%% With theta
display('With theta')

data_path = [data_folder 'Porcine-1_ses-mri01_small_sample_' num2str(nSample) '_concatenate_signal_polyfit_cartesian_all_orientations_with_theta_2_ref.nii.gz'];

% Load data 
data_signal_values = load_nii_img_only(data_path);

temp_shape = size(data_signal_values);
dims = temp_shape(1:3);
data_signal_reshape = permute(data_signal_values,[4 1 2 3]);
data_avg = 0;
for k=1:dims(1)
    for l=1:dims(2)
        for m=1:dims(3)
            if mask_erode(k, l, m) == 1
                subplot(223)
                plot(data_signal_reshape(:,k,l,m))
                
                data_avg = data_avg + data_signal_reshape(:,k,l,m);
                hold on
            end
        end
    end
end

data_avg = data_avg / sum(mask_erode, 'all');
subplot(224)
plot(data_avg)
hold on
% load dictionary with theta
tic()
dict_name = 'SignalWithNoise0_1rep_9orientations_Porcine1_fix_xa_polyfit_cartesian_with_theta.h5py';
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

data_signal_selection = data_avg;

data_signal_replicate = repmat(data_signal_selection,1,dict_length);
mse = sum((data_signal_replicate - dict_signal_reshape).^2, 1);

[mse_max, mse_sort] = sort(mse);

% figure
% plot(data_signal_selection, 'k')
% hold on

nb_map = 1;
for nb = 1:nb_map
    best_dict_match = dict_signal_reshape(:, mse_sort(nb));
    
    FVF(nb) = FVF_reshape(mse_sort(nb));
    gRatio(nb) = gRatio_reshape(mse_sort(nb));
    xi(nb) = xi_reshape(mse_sort(nb));
    T2out(nb) = T2out_reshape(mse_sort(nb));
    T2myel(nb) = T2myel_reshape(mse_sort(nb));
    weight(nb) = weight_reshape(mse_sort(nb));
    
    plot(best_dict_match)
    legend('data', 'dict best match')

    title(['FVF : ' num2str(FVF(nb)) ', gRatio : ' num2str(gRatio(nb)) ', xi : ' num2str(xi(nb)) ...
         ', T2out : ' num2str(T2out(nb)) ', T2myel : ' num2str(T2myel(nb)) ', weight : ' num2str(weight(nb))])
    
end
end

