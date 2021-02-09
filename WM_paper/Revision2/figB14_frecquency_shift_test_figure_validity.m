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

common_model_parameters.phi = 0;

% Classic model
classic_model_parameters = common_model_parameters;

classic_model_parameters.myelin.T2 = 15*1e-3; %s
classic_model_parameters.myelin.weight = 0.5; 
classic_model_parameters.myelin.xi = -0.1; % ppm
classic_model_parameters.myelin.xa = -0.1; % ppm

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

N = 400;
it = 0;

% for l = [0.2 0.3 0.4 0.5]
theta_list = linspace(0, pi/2, 10);
% theta_list = 0;

% theta_list = [0, 2*pi/9, pi/2];
% theta_list = [2*pi/9];
% theta_list = linspace(0, pi/2, 2);
% theta_list = [0 pi/2];
% theta_list = 0;

proportion_list = 0.7;

nb_layer = 10;
myelin_water_proportion = 0.5;
[classic_model, new_model, mask] = axonWithMyelinWater(FVF, gRatio, nb_layer, myelin_water_proportion, N);
classic_model_parameters.dims = size(classic_model);

%%%%% Compute classic model
for kPorportion = 1:length(proportion_list)
    proportion = proportion_list(kPorportion);
    display(['water proportion: ' num2str(proportion)]);
    
    %%%% Models parameter related to myelin water volume fraction
    classic_model_parameters.myelin.weight = proportion;
    myelin_water_proportion = classic_model_parameters.myelin.weight;
    
    %%%%% Original model tensor
    [original_tensor_X, original_total_X] = create2DTensorXFromOneAxon(classic_model, classic_model_parameters);
           
    for kTheta = 1:length(theta_list)
        %%%% Parameter related to theta
        tic()
        theta = theta_list(kTheta);
        common_model_parameters.theta = theta;
        classic_model_parameters.theta = theta;
        new_model_parameters.theta = theta;
        
        display(['theta: ' num2str(theta)]);
        
        common_model_parameters.field_direction = [sin(common_model_parameters.theta)*cos(common_model_parameters.phi), sin(common_model_parameters.theta)*sin(common_model_parameters.phi), cos(common_model_parameters.theta)];
        classic_model_parameters.current_dir = common_model_parameters.field_direction;
        new_model_parameters.field_direction = common_model_parameters.field_direction;
        
        %%%%% Original model field
        original_field_complex = createFieldFrom2DTensorX(original_tensor_X, common_model_parameters);
        original_field = real(original_field_complex); 

        figure
        imagesc(original_field)
        title('original field')
        colorbar
  
        %%%%% Lorentzian correction
        susceptibility_Z = createSusceptibilityZ(original_total_X, classic_model, classic_model_parameters);
        corrected_field = applyFieldLorentzianCavityCorrection(original_field, susceptibility_Z, classic_model_parameters);

        classic_model_parameters.dims = size(classic_model);
        
        options.edges = -10 : 0.25 : 30;
        options.create_figure = 0;
        %             options.fontSize = 20;
        options.LineWidth = 2;
        options.plot= 0;
        
        [~, original_mean_shift(kPorportion, kTheta), original_pos_shift(kPorportion, kTheta),  original_neg_shift(kPorportion, kTheta)] = createHistogramFieldPerturbation(classic_model, original_field, options);
  
%         figure
%         imagesc(corrected_field)
%         title('corrected field')
%         colorbar
        [~, corrected_mean_shift(kPorportion, kTheta), corrected_pos_shift(kPorportion, kTheta), corrected_neg_shift(kPorportion, kTheta)] = createHistogramFieldPerturbation(classic_model, corrected_field, options);
        
        toc()
    end
end


%%%% Compute model with myelin layers
for kPorportion = 1:length(proportion_list)
    proportion = proportion_list(kPorportion);
    display(['water proportion: ' num2str(proportion)]);
    
    %%%% Models parameter related to myelin water volume fraction
    myelin_water_proportion = classic_model_parameters.myelin.weight;
    new_model_parameters.myelin_phospholipid.xi = -0.1 / (1 - proportion);  % myelin anisotropic susceptibility (ppm)
    new_model_parameters.myelin_phospholipid.xa = -0.1 / (1 - proportion); % myelin anisotropic susceptibility (ppm)    
    
    it = it+1;
    display(['number of layers: ' num2str(nb_layer)]);
    tic()
    
    %%%%% Model with layers
    [classic_model, new_model, mask] = axonWithMyelinWater(FVF, gRatio, nb_layer, myelin_water_proportion, N);
    
    for kTheta = 1:length(theta_list)
        %%%% Parameter related to theta
        theta = theta_list(kTheta);
        common_model_parameters.theta = theta;
        new_model_parameters.theta = theta;
        
        display(['theta: ' num2str(theta)]);
        
        common_model_parameters.field_direction = [sin(common_model_parameters.theta)*cos(common_model_parameters.phi), sin(common_model_parameters.theta)*sin(common_model_parameters.phi), cos(common_model_parameters.theta)];
        new_model_parameters.field_direction = common_model_parameters.field_direction;
        
        new_model_parameters.mask = mask;
        
        %%%%% Myelin water layer
        myelin_layer_tensor_X = create2DTensorXFromOneAxonWithMyelinWater(new_model, new_model_parameters);
        
        myelin_layer_field_complex = createFieldFrom2DTensorX(myelin_layer_tensor_X, common_model_parameters);
        myelin_layer_field = real(myelin_layer_field_complex);
        
        figure
        imagesc(myelin_layer_field)
        title('myelin layer field')
        colorbar
                
        options.edges = -10 : 0.25 : 30;
        options.LineWidth = 2;
        options.plot = 1;
        [~, myelin_layer_mean_shift(kPorportion, kTheta), myelin_layer_pos_shift(kPorportion, kTheta), myelin_layer_neg_shift(kPorportion, kTheta)] = createHistogramFieldPerturbationWithMyelinWater(new_model, myelin_layer_field, options);
        
        toc()
    end
end


for kPorportion = 1:length(proportion_list)
    
    x = round(theta_list*180/pi);
    
    figure
    hold on
    plot(x, original_mean_shift(kPorportion, :),'LineWidth',2);
    plot(x, corrected_mean_shift(kPorportion, :),'LineWidth',2);
    plot(x, myelin_layer_mean_shift(kPorportion, :),'LineWidth',2);
           
%     leg = legend('original', 'corrected')
    leg = legend('original', 'corrected', 'myelin layer')
         
    ylabel('mean frequency shift')
    xlabel('\theta')
    
    leg.NumColumns = 2;
    
    set(gca, 'FontSize', 14, 'FontWeight','bold' )
    
end





