clear
% close all

dict_folder = '/project/3015069.04/dictionaries/multi_orientations/BrainSample2/';

data_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals/';
parameter_folder = [data_folder 'parameter_maps/SignalWithNoise05_8rep_with_theta_9orientations_BrainSample2_fix_xa_with_normalization/fa-05/'];

cd(data_folder)

dict_name = 'SignalWithNoise0_1rep_with_theta_9orientations_BrainSample2_fix_xa.h5py';
dict_path = [dict_folder dict_name];
h5disp(dict_path);

FVF_grid =    h5read(dict_path, '/FVFValues');
gRatio_grid = h5read(dict_path, '/gRatioValues');
xi_grid =     h5read(dict_path, '/xiValues');
T2myel_grid = h5read(dict_path, '/T2myelValues');
T2out_grid =  h5read(dict_path, '/T2outValues');
weight_grid =  h5read(dict_path, '/weightValues');
directionsValues =  h5read(dict_path, '/directionsValues');

Signal_values = h5read(dict_path, '/SignalValues');

FVF_grid = squeeze(FVF_grid(:,:,:,1,:,:,:));
gRatio_grid = squeeze(gRatio_grid(:,:,:,1,:,:,:));
xi_grid = squeeze(xi_grid(:,:,:,1,:,:,:));
T2myel_grid = squeeze(T2myel_grid(:,:,:,1,:,:,:));
T2out_grid = squeeze(T2out_grid(:,:,:,1,:,:,:));
weight_grid = squeeze(weight_grid(:,:,:,1,:,:,:));

directions_on_sphere =  squeeze(directionsValues(:,1,1,1,:,1,1,1))';

Signal_values = squeeze(Signal_values);

gRatio_grid = round(100*gRatio_grid)/100; % To obtain a grid 
 
tic()
% Load brain parameter 
prefix = 'BrainSample-2_ses-03_';
suffix = '_fa-05_fix_xa';

FVF_nii = load_untouch_nii([parameter_folder prefix 'FVF' suffix '.nii.gz']);
FVF_brain = FVF_nii.img;
dims = size(FVF_brain);
gRatio_brain = load_nii_img_only([parameter_folder prefix 'gRatio' suffix '.nii.gz']);
xi_brain = load_nii_img_only([parameter_folder prefix 'xi' suffix '.nii.gz']);
T2out_brain = load_nii_img_only([parameter_folder prefix 'T2out' suffix '.nii.gz']);
T2myel_brain = load_nii_img_only([parameter_folder prefix 'T2myel' suffix '.nii.gz']);
weight_brain = load_nii_img_only([parameter_folder prefix 'weight' suffix '.nii.gz']);

mask_path = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_gre_orientation-4_brain_mask.nii.gz';
mask = load_nii_img_only(mask_path);

mask_indices = find(mask == 1);

figure(2)
FVF_list = FVF_brain(mask_indices);
FVF_list(FVF_list < 0.401) = 0.401;
FVF_list(FVF_list > 0.8) = 0.8;
subplot(231)
histogram(FVF_list)
title('FVF')

gRatio_list = gRatio_brain(mask_indices);
gRatio_list(gRatio_list < 0.501) = 0.501;
gRatio_list(gRatio_list > 0.849) = 0.849;
subplot(232)
histogram(gRatio_list)
title('gRatio')

xi_list = xi_brain(mask_indices);
xi_list(xi_list < -0.199) = -0.199;
xi_list(xi_list >  0.199) =  0.199;
subplot(233)
histogram(xi_list)
title('xi')

T2out_list = T2out_brain(mask_indices);
T2out_list(T2out_list < 0.0201) = 0.0201;
T2out_list(T2out_list > 0.0799) = 0.0799;
subplot(235)
histogram(T2out_list)
title('T2out')

T2myel_list = T2myel_brain(mask_indices);
T2myel_list(T2myel_list < 0.00401) = 0.00401;
T2myel_list(T2myel_list > 0.0199) = 0.0199;
subplot(236)
histogram(T2myel_list)
title('T2myel')

weight_list = weight_brain(mask_indices);
weight_list(weight_list < 0.501) = 0.501;
weight_list(weight_list > 4.99) = 4.99;
subplot(236)
histogram(weight_list)
title('weight')

% Load brain dti
dti_path = '/project/3015069.01/derived/BrainSample-2/ses-03/dwi/orientation-9/dti/BrainSample-2_ses-03_dwi_orientation-9_NSA-16_log-euclidean_dti_V1_2_BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-01.nii.gz';
dti = load_nii_img_only(dti_path);
% 
dti_reshape = reshape(dti, [prod(dims) 3]);
dti_list = squeeze(dti_reshape(mask_indices, :));

nb_pixel = length(mask_indices);

length_one_rotation = 23;
nb_of_data_orientations = 9;
lvec = length_one_rotation * nb_of_data_orientations;
total = zeros(nb_pixel,lvec);

signal_first_interpolation = zeros(nb_pixel, lvec, 20);

for k = 1:lvec
    k
    tic()

    for l = 1 : 20
        signal_first_interpolation(:,k,l) = interpn(FVF_grid, gRatio_grid, xi_grid, ...
            T2myel_grid, T2out_grid, weight_grid, ...
            squeeze(Signal_values(k,:,:,:,l,:,:,:)), ...
            FVF_list, gRatio_list, xi_list, T2myel_list, T2out_list, weight_list, ...
            'linear');
    end
    toc()

end
% 
FVF_nii.hdr.dime.dim(1) = 4;
FVF_nii.hdr.dime.dim(5) = lvec;


signal_interpolated_list = zeros([prod(dims) lvec]);
nb_orientations = 3;

for p = 1:nb_pixel
    p
    current_orientation = dti_list(p, :);
    [weights, closest_orientations] = interpolateOrientation(current_orientation, directions_on_sphere, nb_orientations);

    if (length(weights) ~= 1)
        signal_interpolated_list(mask_indices(p), :) = sum(weights.*squeeze(signal_first_interpolation(p,:,closest_orientations)),2);
    else
        signal_interpolated_list(mask_indices(p), :) = signal_first_interpolation(p,:,closest_orientations);
    end
    if ((signal_interpolated_list(mask_indices(p), 1) == 0) || isnan(signal_interpolated_list(mask_indices(p),1)))
        error('Nan values')
    end

end

signal_interpolated_reshape = reshape(signal_interpolated_list, [dims lvec]);

% for j = 1:lvec
%     tic()
%     signal_interpolated_smoothed(:,:,:,j) = smooth3(signal_interpolated_reshape(:,:,:,j), 'gaussian', [5 5 5]);
%     toc()
% end

FVF_nii.img = signal_interpolated_reshape;

FVF_nii.hdr.dime.dim(1) = 4;
FVF_nii.hdr.dime.dim(5) = lvec;

save_untouch_nii(FVF_nii, [parameter_folder 'signal_recovered_noise05_fa-05.nii.gz'])