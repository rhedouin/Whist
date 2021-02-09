clear
% close all
base_folder = '/project/3015069.04/data/3DEM/disp004/';
cd(base_folder)
% For one theta
theta_list = [0, 15, 30, 45, 60, 75, 90];
% theta_list = [15, 30, 45, 60, 75, 90];

% load([base_folder 'tensor_X_3D_disp004.mat'])
load([base_folder 'Model_3D_disp004.mat'])
% keyboard;
load([base_folder 'Mask_3D_disp004.mat'])
load([base_folder 'Mask_3D_disp004_v2.mat'])

for k = 1:length(theta_list)
    theta_degree = theta_list(k);
    display(['theta_degree: ' num2str(theta_degree)])

    theta = deg2rad(theta_degree);
    phi = 0;
    load(['B_' num2str(theta_degree) '_adj_disp004.mat'])
  
    figure(k)
    options.keep_figure = 1;
    options.edges = (-15:0.2:15);
    options.mask = Mask_3D;
    hist_3D{k} = createHistogramFieldPerturbation(Model_3D, B_adj, options);
  
    options.mask = Mask_3D_v2;
    hist_3D{k} = createHistogramFieldPerturbation(Model_3D, B_adj, options);
  
%     save(['B_' num2str(theta_degree) '_adj_disp004'], 'B_adj')
    toc()
    clear B_adj
end
