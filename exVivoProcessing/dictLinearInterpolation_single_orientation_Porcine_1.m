clear
% close all

N = 400;
dico_folder = ['/project/3015069.04/dictionary/'];
signal_folder = [dico_folder 'signals/'];
cd(signal_folder)

signal_name = 'SignalWithNoise0_1rep_fix_rho_with_theta_1orientation.h5py';
signal_path = [signal_folder signal_name];
h5disp(signal_path);

Signal_values = h5read(signal_path, '/SignalValues');

FVF_grid =    h5read(signal_path, '/FVFValues');
gRatio_grid = h5read(signal_path, '/gRatioValues');
xi_grid =     h5read(signal_path, '/xiValues');
theta_grid =  h5read(signal_path, '/thetaValues');
T2myel_grid = h5read(signal_path, '/T2myelValues');
T2out_grid =  h5read(signal_path, '/T2outValues');

[~, theta_order] = sort(squeeze(theta_grid(1,1,1,:,1,1)));

Signal_values = Signal_values(:,:,:,:,theta_order,:,:);

FVF_grid = FVF_grid(:,:,:,theta_order,:,:);
gRatio_grid = gRatio_grid(:,:,:,theta_order,:,:);
xi_grid = xi_grid(:,:,:,theta_order,:,:);
theta_grid = theta_grid(:,:,:,theta_order,:,:);
T2myel_grid = T2myel_grid(:,:,:,theta_order,:,:);
T2out_grid = T2out_grid(:,:,:,theta_order,:,:);

% FVF_grid = squeeze(FVF_grid(:,:,:,1,:,:));
% gRatio_grid = squeeze(gRatio_grid(:,:,:,1,:,:));
% xi_grid = squeeze(xi_grid(:,:,:,1,:,:));
% T2myel_grid = squeeze(T2myel_grid(:,:,:,1,:,:));
% T2out_grid = squeeze(T2out_grid(:,:,:,1,:,:));

% directions_on_sphere =  squeeze(directionsValues(:,1,1,1,:,1,1))';

% Signal_values = squeeze(Signal_values);

gRatio_grid = round(100*gRatio_grid)/100; % To obtain a grid 
 
base_folder = '/project/3015069.01/model/data/BrainSample2/ses-03/';
gre_folder = [base_folder 'gre/'];

for kori = 1:1
    orientation = ['orientation-' num2str(kori)];
    parameter_folder = [gre_folder 'parameter_maps/single_orientation_with_input_theta_with_output_theta/' orientation '_output_theta/'];

tic()
% Load brain parameter 
suffix = '';
FVF_nii = load_untouch_nii([parameter_folder 'FVF' suffix '.nii.gz']);
FVF_brain = FVF_nii.img;
dims = size(FVF_brain);

gRatio_brain = load_nii_img_only([parameter_folder 'gRatio' suffix '.nii.gz']);
xi_brain = load_nii_img_only([parameter_folder 'xi' suffix '.nii.gz']);
T2out_brain = load_nii_img_only([parameter_folder 'T2out' suffix '.nii.gz']);
T2myel_brain = load_nii_img_only([parameter_folder 'T2myel' suffix '.nii.gz']);
theta_brain = load_nii_img_only([parameter_folder 'theta' suffix '.nii.gz']);

mask = load_nii_img_only([gre_folder 'mask_intersection.nii.gz']);
mask_indices = find(mask == 1);

max_FVF = max(FVF_grid(:));
min_FVF = min(FVF_grid(:));

max_gRatio = max(gRatio_grid(:));
min_gRatio = min(gRatio_grid(:));

max_xi = max(xi_grid(:));
min_xi = min(xi_grid(:));

max_T2out = max(T2out_grid(:));
min_T2out = min(T2out_grid(:));

max_T2myel = max(T2myel_grid(:));
min_T2myel = min(T2myel_grid(:));

max_theta = max(theta_grid(:));
min_theta = min(theta_grid(:));

figure(2)
FVF_list = FVF_brain(mask_indices);
FVF_list(FVF_list < min_FVF + eps) = min_FVF + eps;
FVF_list(FVF_list > max_FVF - eps) = max_FVF - eps;
subplot(231)
histogram(FVF_list)
title('FVF')

gRatio_list = gRatio_brain(mask_indices);
gRatio_list(gRatio_list < min_gRatio - eps) = min_gRatio - eps;
gRatio_list(gRatio_list > max_gRatio - eps) = max_gRatio - eps;
subplot(232)
histogram(gRatio_list)
title('gRatio')

xi_list = xi_brain(mask_indices);
xi_list(xi_list < min_xi + eps) = min_xi + eps;
xi_list(xi_list > max_xi - eps) =  max_xi - eps;
subplot(233)
histogram(xi_list)
title('xi')

T2out_list = T2out_brain(mask_indices);
T2out_list(T2out_list < min_T2out + eps) = min_T2out + eps;
T2out_list(T2out_list > max_T2out - eps) = max_T2out - eps;
subplot(234)
histogram(T2out_list)
title('T2out')

T2myel_list = T2myel_brain(mask_indices);
T2myel_list(T2myel_list < min_T2myel + eps) = min_T2myel + eps;
T2myel_list(T2myel_list > max_T2myel - eps) = max_T2myel - eps;
subplot(235)
histogram(T2myel_list)
title('T2myel')

theta_list = theta_brain(mask_indices);
theta_list(theta_list < min_theta + eps) = min_theta + eps;
theta_list(theta_list > max_theta - eps) = max_theta - eps;
subplot(236)
histogram(theta_list)
title('theta')


% % Load brain dti
% dti_file = 'BrainSample-2_ses-03_dwi_orientation-9_NSA-16_log-euclidean_dti_V1_2_BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-01.nii.gz';
% dti_path = [base_folder 'dti/'  dti_file];
% dti = load_nii_img_only(dti_path);
% % 
% dti_reshape = reshape(dti, [prod(dims) 3]);
% dti_list = squeeze(dti_reshape(mask_indices, :));

nb_pixel = length(mask_indices);

length_one_rotation = 23;
nb_of_data_orientations = 1;
lvec = length_one_rotation * nb_of_data_orientations;
total = zeros(nb_pixel,lvec);

% signal_first_interpolation = zeros(nb_pixel, lvec);
signal_list_interpolation = zeros(prod(dims), lvec);

for j = 1:lvec
    j
    tic()

    signal_list_interpolation(mask_indices,j) = interpn(FVF_grid, gRatio_grid, xi_grid, ...
            theta_grid, T2myel_grid, T2out_grid, ...
            squeeze(Signal_values(j,:,:,:,:,:,:)), ...
            FVF_list, gRatio_list, xi_list, theta_list, T2myel_list, T2out_list, ...
            'linear');
    toc()

end
% 

signal_interpolated_reshape = reshape(signal_list_interpolation, [dims lvec]);

FVF_nii.img = signal_interpolated_reshape;
FVF_nii.hdr.dime.dim(1) = 4;
FVF_nii.hdr.dime.dim(5) = lvec;

save_untouch_nii(FVF_nii, [parameter_folder 'signal_recovered_' orientation '.nii.gz'])

end
