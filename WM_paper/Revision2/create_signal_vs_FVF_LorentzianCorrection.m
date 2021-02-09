% This example simulates field perturbation and corresponding multi GRE
% signal from WM models
% Two single models are provided, one 2D model (by defaut) and one 3D model
% (that you can uncomment to run)
% You can load your own WM model create by createOneWMModelExample.m
 
clear
close all


% myelin (required: T2, xi, xa)
model_parameters.myelin.T2 = 15*1e-3;
% model_parameters.myelin.T1 = 500*1e-3;
model_parameters.myelin.proton_density= 0.5; 

model_parameters.myelin.xi = -0.1;  % myelin anisotropic susceptibility (ppm)
model_parameters.myelin.xa = -0.1;  % myelin isotropic susceptibility (ppm)

% intra axonal (required: T2)
model_parameters.intra_axonal.T2 = 50*1e-3;
% model_parameters.intra_axonal.T1 = 1.5;
model_parameters.intra_axonal.proton_density= 1; 
model_parameters.intra_axonal.xi= 0; 

% extra axonal (required: T2)
model_parameters.extra_axonal.T2 = 50*1e-3;
% model_parameters.extra_axonal.T1 = 1.5;
model_parameters.extra_axonal.proton_density= 1; 
model_parameters.extra_axonal.xi= 0; 

% TE (required)
% model_parameters.TE = (2:3:59)*1e-3;
model_parameters.TE = linspace(0.0001,0.08,50); 
nb_TE = length(model_parameters.TE);

% optional, needed to include T1 effect in signal weights
model_parameters.flip_angle = 20;
model_parameters.TR = 60*1e-3;
% 
model_parameters.include_proton_density = 1;
model_parameters.include_T1_effect = 0;

% main magnetic field strength in Tesla (required)
model_parameters.B0 = 3;
% magnetic field orientation (required)

base_model_folder = '/project/3015069.04/WM_Models/N400/';

model_list = [1, 3, 5, 7];
FVF_list = [10, 20, 30, 40, 50, 60, 70, 80];
theta_list = linspace(0, pi/2, 6);
nb_orientations = length(theta_list);

gRatio_round = 60;
fig = figure
colors = linspecer(length(FVF_list));

for kFVF = 1:length(FVF_list)
    FVF_round = FVF_list(kFVF);
    display(['FVF: ' num2str(FVF_round)])
    for num = model_list
        display(['model: ' num2str(num)])

        model_path = [base_model_folder 'FVF' num2str(FVF_round) '_N400_train' num2str(num) '/FVF' num2str(FVF_round) '_gRatio' num2str(gRatio_round) '_N400_train' num2str(num)'];
        load(model_path)
        
        model_parameters.mask = mask;

        number_dims = ndims(mask);
        model_parameters.dims = size(mask);
        
        total_signal = [];
        total_signal_lor = [];
        
        for k = 1:nb_orientations
            
            theta = theta_list(k);
            phi = 0;
            model_parameters.theta = theta;
            model_parameters.phi = phi;
            model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
            model_parameters.current_dir = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
            
%             orientations(k,:) = model_parameters.field_direction ;
            % Estimate the relative weights of each compartment
            model_parameters = computeCompartmentSignalWeight(model_parameters);
                        
            %%%%%%%%% Simulate the field perturbation from the WM model and the multi GRE signals
            [signal_original, field] = simulateSignalFromModel(axon_collection, model_parameters);
            [signal_original_lor, field_lor] = simulateSignalFromModelWithLorentzianCorrection(axon_collection, model_parameters);

            total_signal = [total_signal signal_original.total_normalized];
            total_signal_lor = [total_signal_lor signal_original_lor.total_normalized];
        end

        figure(1)
        subplot(211)
        plot(real(total_signal), '-', 'Color', colors(kFVF, :), 'LineWidth', 1.5);
        ylabel('real')
        hold on 
        
        subplot(212)
        p(kFVF) = plot(imag(total_signal), '-', 'Color', colors(kFVF, :), 'LineWidth', 1.5);
        ylabel('imag')      
        hold on
        
        
        
        figure(2)
        subplot(211)
        plot(real(total_signal_lor), '-', 'Color', colors(kFVF, :), 'LineWidth', 1.5);
        ylabel('real')
        hold on 
        
        subplot(212)
        p(kFVF) = plot(imag(total_signal_lor), '-', 'Color', colors(kFVF, :), 'LineWidth', 1.5);
        ylabel('imag')      
        hold on
    end    
end

figure(2)
leg = legend(p, '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8');
title(leg, 'FVF');
leg.Location = 'southwest';
leg.NumColumns = 2;

subplot(212)
for k = 1:nb_orientations -1
    vline(nb_TE * k, '--k')
end


subplot(211)
for k = 1:nb_orientations -1
    vline(nb_TE * k, '--k')
end

subplot(211)
set(gca, 'FontSize', 20, 'FontWeight','bold', 'XTick', [])

subplot(212)
xlabel('concatenate signal')
set(gca, 'FontSize', 20, 'FontWeight','bold', 'XTick', [])








