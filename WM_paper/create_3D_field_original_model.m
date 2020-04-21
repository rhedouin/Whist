clear
close all
base_folder = '/project/3015069.04/data/3DEM/Dispersion004/';
cd(base_folder)
% For one theta
theta_list = [0, 15, 30, 45, 60, 75, 90];
theta_list = [15, 30, 45, 60, 75, 90];

load([base_folder 'tensor_X_3D_original.mat'])
load([base_folder 'Mask_3D_original.mat'])

for k = 1:length(theta_list)
    theta_degree = theta_list(k);
    display(['theta_degree: ' num2str(theta_degree)])

    theta = deg2rad(theta_degree);
    phi = 0;
    model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    model_parameters.B0 = 3;
    tic()
    B_adj = createFieldFrom3DTensorX(tensor_X_3D, model_parameters);
  
    B_adj = real(B_adj);    
    B_adj = single(B_adj.*Mask_3D);
    toc()
    tic()
    save(['B_' num2str(theta_degree) '_adj_disp004'], 'B_adj')
    toc()
    clear B_adj
end
