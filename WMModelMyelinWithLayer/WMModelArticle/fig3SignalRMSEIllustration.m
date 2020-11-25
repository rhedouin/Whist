% test Axon Function 
close all;
clear;

%%%%%%%%% Set parameters
% Common to classic and new
% intra axonal (required: T2)
common_model_parameters.intra_axonal.T2 = 50*1e-3;
common_model_parameters.intra_axonal.weight = 1;
common_model_parameters.intra_axonal.xi= 0; 

% extra axonal (required: T2)
common_model_parameters.extra_axonal.T2 = 50*1e-3;
common_model_parameters.extra_axonal.weight = 1;
common_model_parameters.extra_axonal.xi= 0; 

common_model_parameters.TE = (1:1:60)*1e-3;
common_model_parameters.field_direction = [1, 0, 0];

common_model_parameters.B0 = 3;

% Classic model
classic_model_parameters = common_model_parameters;

classic_model_parameters.myelin.T2 = 15*1e-3;
classic_model_parameters.myelin.weight = 0.5; 
classic_model_parameters.myelin.xi = -0.1; 
classic_model_parameters.myelin.xa = -0.1; 

% New model
new_model_parameters = common_model_parameters;

new_model_parameters.myelin_phospholipid.T2 = 0.1*1e-3;
new_model_parameters.myelin_phospholipid.weight= 1; 
new_model_parameters.myelin_phospholipid.xi = -0.1 / (1 - classic_model_parameters.myelin.weight);  % myelin anisotropic susceptibility (ppm)
new_model_parameters.myelin_phospholipid.xa = -0.1 / (1 - classic_model_parameters.myelin.weight);  % myelin anisotropic susceptibility (ppm)

new_model_parameters.myelin_water.T2 = 15*1e-3;
new_model_parameters.myelin_water.weight= 1; 
new_model_parameters.myelin_water.xi = 0; 

FVF = 0.4;
gRatio = 0.6;

N = 1000;
it = 0;

% for l = [0.2 0.3 0.4 0.5]
theta_list = linspace(0, pi/2, 19);
proportion_list = [0.2 0.3 0.4 0.5];
nb_layer_list = [0 5 10 15 20];

phi = 0;

for kPorportion = 1:length(proportion_list)
    proportion = proportion_list(kPorportion);
    display(['water proportion: ' num2str(proportion)]);

    for kTheta = 1:length(theta_list)
        theta = theta_list(kTheta);
        display(['theta: ' num2str(theta)]);

        common_model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
        classic_model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
        new_model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
         
        classic_model_parameters.myelin.weight = proportion;
        myelin_water_proportion = classic_model_parameters.myelin.weight;
        new_model_parameters.myelin_phospholipid.xi = -0.1 / (1 - proportion);  % myelin anisotropic susceptibility (ppm)
        new_model_parameters.myelin_phospholipid.xa = -0.1 / (1 - proportion); % myelin anisotropic susceptibility (ppm)
        
        for kLayer = 1:length(nb_layer_list)
            nb_layer = nb_layer_list(kLayer);
            
            it = it+1;
            display(['number of layers: ' num2str(nb_layer)]);
            tic()
            
            [classic_model, new_models, mask] = axonWithMyelinWater(FVF, gRatio, nb_layer, myelin_water_proportion, N);

            new_model_parameters.mask = mask;
            classic_model_parameters.mask = mask;
            
            classic_model_parameters.dims = size(classic_model);
            
            if nb_layer == 0
                tensor_X = create2DTensorXFromOneAxon(classic_model, classic_model_parameters);
            else
                tensor_X = create2DTensorXFromOneAxonWithMyelinWater(new_models, new_model_parameters);
            end
            
            field_complex = createFieldFrom2DTensorX(tensor_X, common_model_parameters);
            field = real(field_complex);
                        
            options.edges = -10 : 0.25 : 60;
            options.fontSize = 20;
            options.LineWidth = 2;
            
            if nb_layer == 0
                signals{kPorportion, kTheta, kLayer} = simulateSignalFromField(classic_model, field, classic_model_parameters);
%                 createHistogramFieldPerturbation(classic_model, field{it}, options);
            else
                signals{kPorportion, kTheta, kLayer} = simulateSignalFromFieldWithMyelinWater(new_models, field, new_model_parameters);
%                 createHistogramFieldPerturbationWithMyelinWater(new_models{it}, field{it}, options);
            end       

            simi_signal(kPorportion, kTheta, kLayer) = computeSignalSimilarityIndex(signals{kPorportion, kTheta, kLayer}.total_normalized, signals{kPorportion, kTheta, 1}.total_normalized)
       
            toc()
        end
    end
    figure
    hold on
    
    x = round(theta_list*180/pi);

    for kLayer = 1:length(nb_layer_list)
        plot(x, simi_signal(kPorportion, :, kLayer))
    end
    legend('0', '5', '10', '15', '20')
end

save('simi_signal_N1000.mat','simi_signal')



