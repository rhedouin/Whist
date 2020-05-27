clear
close all

base_folder = '/project/3015069.04/data/';
twoD_folder = [base_folder  '2DRM/'];

ext = 'disp04';
threeD_folder = [base_folder  '3DEM/' ext '/'];

cd(base_folder)

% 2D creation
%%%%%%%%%%% Set parameters

% myelin (required: T2, xi, xa)
model_parameters.myelin.T2 = 16*1e-3;
model_parameters.myelin.proton_density= 0.5; 

model_parameters.myelin.xi = -0.1;  % myelin anisotropic susceptibility (ppm)
model_parameters.myelin.xa = -0.1;  % myelin isotropic susceptibility (ppm)

% intra axonal (required: T2)
model_parameters.intra_axonal.T2 = 60*1e-3;
model_parameters.intra_axonal.proton_density= 1; 
model_parameters.intra_axonal.xi= 0; 

% extra axonal (required: T2)
model_parameters.extra_axonal.T2 = 60*1e-3;
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
model_parameters.myelin.weight = 2;

model_parameters.nb_orientation_for_dispersion = 100;

% model_parameters.kappa_list = [10000 18 9 5.5 3.5 2 1];
% model_parameters.dispersion_list = [0.001 0.1 0.2 0.3 0.4 0.5 0.6];

model_parameters.kappa_list = [10000 9 3.5 1];
model_parameters.dispersion_list = [0.001 0.2 0.4 0.6];

% model_parameters.kappa_list = [10000 5.5 1];
% model_parameters.dispersion_list = [0.001 0.3 0.6];
% 
% model_parameters.kappa_list = [1];
% model_parameters.dispersion_list = [0.6];

% For one theta
theta_list = [0, 15, 30, 45, 60, 75, 90];
% theta_list = [0];

% theta_list = [45];

nb_models = 10;

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
        Mask_3D = single(Mask_3D);
        
    AVF = length(find(Model_3D.*Mask_3D == 0.5));
    MVF = length(find(Model_3D.*Mask_3D == 1));
    
    FVF = (AVF + MVF) / sum(Mask_3D, 'all');
    g_ratio = sqrt(AVF / (AVF + MVF));
            
    load([threeD_folder 'B_' num2str(theta_degree) '_adj_' ext '.mat'])
        
%     figure(k)
    options.keep_figure = 0;
    options.edges = (-15:0.2:15);
    options.mask = Mask_3D;
    hist_3D{k} = createHistogramFieldPerturbation(Model_3D, B_adj, options);
    
    model_parameters.mask = Mask_3D;
    
    signal_3D{k} =  simulateSignalFromField(Model_3D, B_adj, model_parameters);

%     figure(k + 10)
%     subplot(211)
%     hold on
%     plot(abs(signal_3D.total_normalized), '-')
%     subplot(212)
%     hold on
%     plot(phase(signal_3D.total_normalized), '-')
    
    clear options Model_3D Mask_3D B_adj
    
    for l = 1:length(model_parameters.dispersion_list)
        
        display('2D model: load and create histogram')
        model_parameters.kappa = model_parameters.kappa_list(l);
        model_parameters.dispersion = model_parameters.dispersion_list(l);
        
        display(['dispersion: ' num2str(model_parameters.dispersion)])
        
        %For one model
        for m = 1:nb_models
            load([twoD_folder '2DModel_FVF054_gRatio_067_v' num2str(m) '.mat'])
            model_2D = model;
            mask_2D = mask;
            clear model mask
            
            model_parameters.mask = mask_2D;
            model_parameters.dims = size(mask_2D);
            
            model_2D_replicate = repmat(model_2D, [1 1 model_parameters.nb_orientation_for_dispersion]);
            mask_2D_replicate = repmat(mask_2D, [1 1 model_parameters.nb_orientation_for_dispersion]);
            
            %%%%%%%%% Simulate the field perturbation from the WM model and the multi GRE signals
            [signal_2D{k,l,m}, field_2D_dispersion] = simulateSignalWithDispersionFromModel(axon_collection, model_parameters);
            
%             if l == 1
%                 figure(k + 10)
%                 subplot(211)
%                 plot(abs(signal_2D{k,l,m}.total_normalized), '--')
%                 subplot(212)
%                 plot(phase(signal_2D{k,l,m}.total_normalized), '--')
%             elseif l == 2
%                 figure(k + 10)
%                 subplot(211)
%                 plot(abs(signal_2D{k,l,m}.total_normalized), '-.')
%                 subplot(212)
%                 plot(phase(signal_2D{k,l,m}.total_normalized), '-.')
%             elseif l == 3
%                 figure(k + 10)
%                 subplot(211)
%                 plot(abs(signal_2D{k,l,m}.total_normalized), ':.')
%                 subplot(212)
%                 plot(phase(signal_2D{k,l,m}.total_normalized), ':.')
%             elseif l == 4
%                 figure(k + 10)
%                 subplot(211)
%                 plot(abs(signal_2D{k,l,m}.total_normalized), 'r--')
%                 subplot(212)
%                 plot(phase(signal_2D{k,l,m}.total_normalized), 'r--')
%             end
            
            options.mask = mask_2D_replicate;
            options.edges = (-15:0.2:15);
            options.keep_figure = 1;
            
%             figure(k)
%             
%             if l == 1
%                 options.line_style = '--';
%             elseif l == 2
%                 options.line_style = '-.';
%             elseif l == 3
%                 options.line_style = ':.';
%             elseif l == 4
%                 options.line_style = '--';
%             end
            
            hist_2D{k,l,m} = createHistogramFieldPerturbation(model_2D_replicate, field_2D_dispersion, options);

            diff_histo(k,l,m) = sqrt(sum((hist_2D{k,l,m}.intra_axonal - hist_3D{k}.intra_axonal).^2 + ...
                           (hist_2D{k,l,m}.extra_axonal - hist_3D{k}.extra_axonal).^2 + ...
                           (hist_2D{k,l,m}.myelin - hist_3D{k}.myelin).^2));
                        
            simi_signal(k,l,m) = computeSignalSimilarityIndex(signal_2D{k,l,m}.total_normalized, signal_3D{k}.total_normalized)
        end
    end
    clear model_2D mask_2D options
end

%%%%%%%%%%%%%%%%%%%%%%%
save([twoD_folder 'signal_and_histo_distance_2D_vs_3D_' ext '_new_mask_100_models_paper_parameter_values_disp04.mat'], 'diff_histo', 'simi_signal', 'signal_2D', 'signal_3D')



