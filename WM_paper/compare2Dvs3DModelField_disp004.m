clear
close all

base_folder = '/project/3015069.04/data/';
twoD_folder = [base_folder  '2DRM/'];

ext = 'disp02';
threeD_folder = [base_folder  '3DEM/' ext '/'];

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

for k = 1:length(theta_list)
    theta_degree = theta_list(k);
    display(['theta_degree: ' num2str(theta_degree)])

    theta = deg2rad(theta_degree);
    phi = 0;
    model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    
    %%%%%%% Compute 3D model histogram
    display('3D model: load and create histogram')
    
    load([threeD_folder 'Mask_3D_' ext '.mat'])
    load([threeD_folder 'Model_3D_' ext '.mat'])
    load([threeD_folder 'B_' num2str(theta_degree) '_adj_' ext '.mat'])
%         
%     figure(k)
%     options.keep_figure = 0;
%     options.edges = (-15:0.2:15);
%     options.mask = Mask_3D;
%     hist_3D{k} = createHistogramFieldPerturbation(Model_3D, B_adj, options);
    
%     clear options Model_3D Mask_3D B_adj
    
    for l = 1:length(model_parameters.dispersion_list)
        
        display('2D model: load and create histogram')
        model_parameters.kappa = model_parameters.kappa_list(l);
        model_parameters.dispersion = model_parameters.dispersion_list(l);
        
        display(['dispersion: ' num2str(model_parameters.dispersion)])
        
        %For one model
        for m = 1:nb_models
            load([twoD_folder '2DModel_FVF049_gRatio_064_v' num2str(m) '.mat'])
            model_2D = model;
            mask_2D = mask;
            clear model mask
            
            model_parameters.mask = mask_2D;
            model_parameters.dims = size(mask_2D);
            
            model_2D_replicate = repmat(model_2D, [1 1 model_parameters.nb_orientation_for_dispersion]);
            mask_2D_replicate = repmat(mask_2D, [1 1 model_parameters.nb_orientation_for_dispersion]);
            
            %%%%%%%%%% Simulate the field perturbation from the WM model and the multi GRE signals
            [~, field_2D_dispersion] = simulateSignalWithDispersionFromModel(axon_collection, model_parameters);

%             [~, field_2D] = simulateSignalFromModel(axon_collection, model_parameters);
            
            figure(k+10)
            subplot(221)
            imagesc(field_2D_dispersion(:,:,1));
            colorbar
            caxis([-5 8])
            subplot(222)
            imagesc(field_2D_dispersion(:,:,2));
            colorbar
            caxis([-5 8])
            subplot(223)
            imagesc(field_2D_dispersion(:,:,3));
            colorbar
            caxis([-5 8])
            subplot(224)
            imagesc(field_2D_dispersion(:,:,4));
            colorbar
            caxis([-5 8])

            options.mask = mask_2D_replicate;          
            options.edges = (-15:0.2:15);
            
            figure(k)
            
            options.keep_figure = 0;
            options.edges = (-15:0.2:15);
            options.mask = Mask_3D;
            hist_3D{k} = createHistogramFieldPerturbation(Model_3D, B_adj, options);
            options.keep_figure = 1;
            %             hist_2D{k,l,m} = createHistogramFieldPerturbation(model_2D_replicate, field_2D_dispersion, options);
            options.mask = mask_2D;
            hist_2D{k,l,m} = createHistogramFieldPerturbation(model_2D_replicate, field_2D_dispersion, options);
         
            
            d(k,l,m) = sum((hist_2D{k,l,m}.intra_axonal - hist_3D{k}.intra_axonal).^2 + ...
                           (hist_2D{k,l,m}.extra_axonal - hist_3D{k}.extra_axonal).^2 + ...
                           (hist_2D{k,l,m}.myelin - hist_3D{k}.myelin).^2)
        end
    end
    clear model_2D mask_2D options
end
              return;

% save([twoD_folder 'distance_l2_2D_vs_3D_' ext '.mat'], 'd')




