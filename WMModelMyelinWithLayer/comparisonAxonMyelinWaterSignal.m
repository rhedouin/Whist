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

new_model_parameters.myelin_phospholipid.T2 = 1*1e-3;
new_model_parameters.myelin_phospholipid.weight= 1; 
new_model_parameters.myelin_phospholipid.xi = -0.2;  % myelin anisotropic susceptibility (ppm)
new_model_parameters.myelin_phospholipid.xa = -0.2;  % myelin anisotropic susceptibility (ppm)

new_model_parameters.myelin_water.T2 = 15*1e-3;
new_model_parameters.myelin_water.weight= 1; 
new_model_parameters.myelin_water.xi = 0; 

FVF = 0.3;
gRatio = 0.7;
% myelin_phospholipid_proportion = 1 - classic_model_parameters.myelin.weight;
% myelin_phospholipid_proportion = 0.1;

N = 1500;

it = 0;

for l = [0.5]
      myelin_water_proportion = classic_model_parameters.myelin.weight;

%     for k = [0 1 4 10 25]
      for k = [0 5 10 15 20]
%     for k = [0 1]

%             for k = [0 10]
        
        it = it+1;
        display(['number of layers: ' num2str(k)]);
        nb_myelin_water_layer = k;
        
        [classic_model, new_model, mask] = axonWithMyelinWater(FVF, gRatio, nb_myelin_water_layer, myelin_water_proportion, N);

        
        figure
        imagesc(new_model)
        colorbar
        colormap hot
        
        new_model_parameters.mask = mask;
        classic_model_parameters.mask = mask;
        
        classic_model_parameters.dims = size(classic_model);
        
        if k == 0
            tensor_X = create2DTensorXFromOneAxon(classic_model, classic_model_parameters);
        else    
            tensor_X = create2DTensorXFromOneAxonWithMyelinWater(new_model, new_model_parameters);
        end
%         
%         figure
%         subplot(231)
%         imagesc(tensor_X(:,:,1))
%         colorbar;
%         
%         subplot(232)
%         imagesc(tensor_X(:,:,2))
%         colorbar;
%         
%         subplot(233)
%         imagesc(tensor_X(:,:,3))
%         colorbar;
%         
%         subplot(234)
%         imagesc(tensor_X(:,:,4))
%         colorbar;
% 
%         subplot(235)
%         imagesc(tensor_X(:,:,5))
%         colorbar;
%         
%         subplot(236)
%         imagesc(tensor_X(:,:,6))
%         colorbar;
% 
%         keyboard;

        field_complex = createFieldFrom2DTensorX(tensor_X, common_model_parameters);
        field{it} = real(field_complex);
        
        figure
        imagesc(field{it})
        colorbar
        caxis([-10 10])
        
        if k == 0
            signals{it} = simulateSignalFromField(classic_model, field{it}, classic_model_parameters);
        else
            signals{it} = simulateSignalFromFieldWithMyelinWater(new_model, field{it}, new_model_parameters);            
        end
        
    
        figure(10)
        subplot(211)
        ylabel('abs')
        plot(abs(signals{it}.total_normalized))
        hold on
        subplot(212)
                ylabel('phase')

        plot(angle(signals{it}.total_normalized))
        hold on
        
        figure
        imagesc(field{it})
    end
    
end

legend('0','5','10','15','20')
return;




return;
figure
subplot(231)
imagesc(tensor_X(:,:,1));
colorbar;

subplot(232)
imagesc(tensor_X(:,:,2));
colorbar;

subplot(233)
imagesc(tensor_X(:,:,3));
colorbar;

subplot(234)
imagesc(tensor_X(:,:,4));
colorbar;

subplot(235)
imagesc(tensor_X(:,:,5));
colorbar;

subplot(236)
imagesc(tensor_X(:,:,6));
colorbar;

