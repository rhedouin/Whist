clear
close all

base_folder = '/project/3015069.04/data/';
twoD_folder = [base_folder  '2DRM/'];

cd(base_folder)

% 2D creation
%%%%%%%%%%% Set parameters

% myelin (required: T2, xi, xa)
model_parameters.myelin.T2 = 15*1e-3;
model_parameters.myelin.proton_density= 0.5; 

model_parameters.myelin.xi = -0.1;  % myelin anisotropic susceptibility (ppm)
model_parameters.myelin.xa = -0.1;  % myelin isotropic susceptibility (ppm)

% intra axonal (required: T2)
model_parameters.intra_axonal.T2 = 50*1e-3;
model_parameters.intra_axonal.proton_density= 1; 
model_parameters.intra_axonal.xi= 0; 

% extra axonal (required: T2)
model_parameters.extra_axonal.T2 = 50*1e-3;
model_parameters.extra_axonal.proton_density= 1; 
model_parameters.extra_axonal.xi= 0; 

% main magnetic field strength in Tesla (required)
model_parameters.B0 = 3;
% magnetic field orientation (required)

% TE (required)
% model_parameters.TE = (2:3:59)*1e-3;
model_parameters.TE = linspace(0.0001,0.08,100); 
% 
model_parameters.include_proton_density = 1;
model_parameters.include_T1_effect = 0;
%%%%%%%%%%%

% Estimate the relative weights of each compartment
model_parameters = computeCompartmentSignalWeight(model_parameters);
model_parameters.intra_axonal.weight = 1;
model_parameters.extra_axonal.weight = 1;
model_parameters.myelin.weight = 1.5;

model_parameters.nb_orientation_for_dispersion = 100;

model_parameters.kappa_list = [10000 18 9 5.5 3.5 2 1];
model_parameters.dispersion_list = [0.001 0.1 0.2 0.3 0.4 0.5 0.6];

% For one theta
theta_list = [0, 15, 30, 45, 60, 75, 90];

% theta_list = [0];

nb_models = 1;

ext_list = {'disp004', 'disp02', 'disp04'};

for l = 1:length(ext_list)
    ext = ext_list{l};
    
    threeD_folder = [base_folder  '3DEM/' ext '/'];

    %%%%%%% Compute 3D model histogram
    display('3D model: load and create histogram')
    
    load([threeD_folder 'Mask_3D_' ext '.mat'])
    load([threeD_folder 'Model_3D_' ext '.mat'])
    
    for k = 1:length(theta_list)
        theta_degree = theta_list(k);
        load([threeD_folder 'B_' num2str(theta_degree) '_adj_' ext '.mat'])

        display(['theta_degree: ' num2str(theta_degree)])
        
        theta = deg2rad(theta_degree);
        phi = 0;
        model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
        
      
        figure(k)
        options.keep_figure = 1;
        options.edges = (-15:0.2:15);
        options.mask = Mask_3D;
        if l == 1
            options.line_style = '-';
        elseif l == 2
            options.line_style = '--';
        else     
            options.line_style = ':';
        end
        
        hist_3D{k} = createHistogramFieldPerturbation(Model_3D, B_adj, options);
        
%         clear options Model_3D Mask_3D B_adj
        
    end
end
              return;

% save([twoD_folder 'distance_l2_2D_vs_3D_' ext '.mat'], 'd')




