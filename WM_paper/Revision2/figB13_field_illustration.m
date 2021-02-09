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

theta = pi/2;
phi = 0;

common_model_parameters.theta = theta;
common_model_parameters.field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
common_model_parameters.current_dir = common_model_parameters.field_direction;

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

for l = [0.5]
    classic_model_parameters.myelin.weight = l;
    myelin_water_proportion = classic_model_parameters.myelin.weight;
    new_model_parameters.myelin_phospholipid.xi = -0.1 / (1 - myelin_water_proportion);  % myelin anisotropic susceptibility (ppm)
    new_model_parameters.myelin_phospholipid.xa = -0.1 / (1 - myelin_water_proportion);  % myelin anisotropic susceptibility (ppm)
    
    k = 5;
    
    it = it+1;
    display(['number of layers: ' num2str(k)]);
    nb_myelin_water_layer = k;
    
    [classic_model, new_model, mask] = axonWithMyelinWater(FVF, gRatio, nb_myelin_water_layer, myelin_water_proportion, N);
    classic_model_parameters.dims = size(classic_model);
    
    new_model_parameters.mask = mask;
    classic_model_parameters.mask = mask;
    
    [signal_corrected, field_corrected, field_original] = simulateSignalForOneAxonlWithLorentzianCorrection(classic_model, classic_model_parameters); 
    
    tensor_X_myelin_layers = create2DTensorXFromOneAxonWithMyelinWater(new_model, new_model_parameters);
    
    field_complex_myelin_layers  = createFieldFrom2DTensorX(tensor_X_myelin_layers , common_model_parameters);
    field_myelin_layers  = real(field_complex_myelin_layers);
    
    figure
    imagesc(field_myelin_layers)
    
    set(gca,'XTick',[], 'YTick', [], 'fontsize', 20, 'fontweight', 'bold')
    
    options.create_figure = 1;
    options.line_style = '--';
    options.mask = mask;
    options.edges = (-15:0.1:15);
    options.LineWidth = 2;
    options.xlim = [-15 15];
%     options.ylim = [0 0.04];
    
    createHistogramFieldPerturbation(classic_model, field_original, options);
    
    options.create_figure = 0;
    options.line_style = '-';
    createHistogramFieldPerturbation(classic_model, field_corrected, options);
    set(gca, 'FontSize', 16, 'FontWeight','bold' )
    
    options.create_figure = 1;
    options.line_style = '-';
 
    createHistogramFieldPerturbationWithMyelinWater(new_model, field_myelin_layers, options);
    set(gca, 'FontSize', 16, 'FontWeight','bold' )

end
