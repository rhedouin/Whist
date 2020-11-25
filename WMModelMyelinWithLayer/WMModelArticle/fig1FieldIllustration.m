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
% myelin_phospholipid_proportion = 1 - classic_model_parameters.myelin.weight;
% myelin_phospholipid_proportion = 0.1;

N = 1000;

it = 0;

for l = [0.4]
      classic_model_parameters.myelin.weight = l;
      myelin_water_proportion = classic_model_parameters.myelin.weight;
      new_model_parameters.myelin_phospholipid.xi = -0.1 / l;  % myelin anisotropic susceptibility (ppm)
      new_model_parameters.myelin_phospholipid.xa = -0.1 / l;  % myelin anisotropic susceptibility (ppm)

      for k = [0 5]
        
        it = it+1;
        display(['number of layers: ' num2str(k)]);
        nb_myelin_water_layer = k;
        
        [classic_model, new_models{it}, mask] = axonWithMyelinWater(FVF, gRatio, nb_myelin_water_layer, myelin_water_proportion, N);
        
        new_model_parameters.mask = mask;
        classic_model_parameters.mask = mask;
        
        classic_model_parameters.dims = size(classic_model);
        
        if k == 0
            tensor_X = create2DTensorXFromOneAxon(classic_model, classic_model_parameters);
        else    
            tensor_X = create2DTensorXFromOneAxonWithMyelinWater(new_models{it}, new_model_parameters);
        end

        field_complex = createFieldFrom2DTensorX(tensor_X, common_model_parameters);
        field{it} = real(field_complex);
        
        figure
        imagesc(field{it})
                      
        set(gca,'XTick',[], 'YTick', [], 'fontsize', 20, 'fontweight', 'bold')
        
        options.edges = -20 : 0.25 : 20;
        options.fontSize = 16;
        options.LineWidth = 2;
        
        if k == 0
            signals{it} = simulateSignalFromField(classic_model, field{it}, classic_model_parameters);
            createHistogramFieldPerturbation(classic_model, field{it}, options);
        else
            signals{it} = simulateSignalFromFieldWithMyelinWater(new_models{it}, field{it}, new_model_parameters);            
            createHistogramFieldPerturbationWithMyelinWater(new_models{it}, field{it}, options);
        end 
    end
end
