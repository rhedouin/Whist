% close all
clear

dict_folder = '/project/3015069.04/dictionaries/multi_orientations/Porcine2/fix_xa_large_FVF_20_directions/';
data_folder = '/project/3015069.01/derived/Porcine-2/ses-mri01/concatenate_signals/individual_samples/';
    
for kSample = 1:1
%     figure(kSample)
%     figure

nSample = kSample;

experience_name = 'Porcine-2';


%%%%%%%%%%%%%%%% With theta
display('With theta')

sample_name = [experience_name '_ses-mri01_small_sample_' num2str(nSample) '_concatenate_signal_all_orientations_polyfit_cartesian_with_theta.nii.gz'];
data_path = [data_folder sample_name];

mask_name = [experience_name '_ses-mri01_small_sample_' num2str(nSample) '_mask.nii.gz'];
mask_path = [data_folder mask_name];
% Load data 
data_signal_values = load_nii_img_only(data_path);
mask = load_nii_img_only(mask_path);

temp_shape = size(data_signal_values);
dims = temp_shape(1:3);
data_signal_reshape = permute(data_signal_values,[4 1 2 3]);

plot_signal = 0;
if plot_signal
    data_avg = 0;
    for k=1:dims(1)
        for l=1:dims(2)
            for m=1:dims(3)
                if  mask(k, l, m) == 1
                    subplot(223)
                    plot(data_signal_reshape(:,k,l,m))
                    
                    data_avg = data_avg + data_signal_reshape(:,k,l,m);
                    hold on
                end
            end
        end
    end
    
    data_avg = data_avg / sum(mask, 'all');
    subplot(224)
    plot(data_avg)
    hold on
end

% load dictionary with theta
tic()
dict_name = 'SignalWithNoise0_1rep_6orientations_12TE_Porcine2_fix_xa_large_FVF_20_directions_polyfit_cartesian_with_theta.h5py';
dict_path = [dict_folder dict_name];
keyboard;
dict_signal_values = h5read(dict_path, '/SignalValues');
dict_signal_shape = size(dict_signal_values);
signal_length = dict_signal_shape(1);
dict_shape = dict_signal_shape(2:end);
dict_length = prod(dict_shape);

FVF_grid =    h5read(dict_path, '/FVFValues');
gRatio_grid = h5read(dict_path, '/gRatioValues');
xiMyelin_grid =     h5read(dict_path, '/xiMyelinValues');
% xaMyelin_grid =     h5read(dict_path, '/xaMyelinValues');
T2IntraExtraAxonal_grid =    h5read(dict_path, '/T2IntraExtraAxonalValues');
T2Myelin_grid=    h5read(dict_path, '/T2MyelinValues');
weight_grid =    h5read(dict_path, '/weightValues');
direction_grid = h5read(dict_path, '/directionsValues');

dict_signal_reshape = reshape(dict_signal_values, signal_length, dict_length);
FVF_reshape = FVF_grid(:);
gRatio_reshape = gRatio_grid(:);
xiMyelin_reshape = xiMyelin_grid(:);
T2IntraExtraAxonal_reshape = T2IntraExtraAxonal_grid(:);
T2Myelin_reshape = T2Myelin_grid(:);
weight_reshape = weight_grid(:);
% 
% data_signal_selection = data_avg;
% 
% data_signal_replicate = repmat(data_signal_selection,1,dict_length);
% mse = sum((data_signal_replicate - dict_signal_reshape).^2, 1);
% 
% [mse_max, mse_sort] = sort(mse);
% 
% % figure
% % plot(data_signal_selection, 'k')
% % hold on

nb_map = 1;
for nb = 1:nb_map
    best_dict_match = dict_signal_reshape(:, mse_sort(nb));
    
    FVF(nb) = FVF_reshape(mse_sort(nb));
    gRatio(nb) = gRatio_reshape(mse_sort(nb));
    xiMyelin(nb) = xiMyelin_reshape(mse_sort(nb));
    T2IntraExtraAxonal(nb) = T2IntraExtraAxonal_reshape(mse_sort(nb));
    T2Myelin(nb) = T2Myelin_reshape(mse_sort(nb));
    weight(nb) = weight_reshape(mse_sort(nb));
    
    plot(best_dict_match)
    legend('data', 'dict best match')

    title(['FVF : ' num2str(FVF(nb)) ', gRatio : ' num2str(gRatio(nb)) ', xi : ' num2str(xi(nb)) ...
         ', T2out : ' num2str(T2out(nb)) ', T2myel : ' num2str(T2myel(nb)) ', weight : ' num2str(weight(nb))])
    
end
end

