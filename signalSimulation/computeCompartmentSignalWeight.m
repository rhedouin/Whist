function model_parameters = computeCompartmentSignalWeight(model_parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function set the weight of each compartment which reflects the water
% signal received.
% It can comports the water proton density and/or the T1 effect. 
% If is advisible to at least include the waterproton density. 

%%%%%%%%%%%%%%%% Input required
% % If include_proton_density
% proton density of each compartment
%
% % If include_T1_effect%  
% T1 of each compartment
% TR
% flip angle 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set defaut parameters
if  ~isfield(model_parameters, 'include_proton_density')
    model_parameters.include_proton_density = 0;
end

if  ~isfield(model_parameters, 'include_T1_effect')
    model_parameters.include_T1_effect = 0;
end

if (~model_parameters.include_proton_density && ~model_parameters.include_T1_effect)
    model_parameters.myelin.weight = 1;
    model_parameters.intra_axonal.weight = 1;
    model_parameters.extra_axonal.weight = 1;
    
elseif (model_parameters.include_proton_density && ~model_parameters.include_T1_effect)
    model_parameters.myelin.weight = model_parameters.myelin.proton_density;
    model_parameters.intra_axonal.weight = model_parameters.intra_axonal.proton_density;
    model_parameters.extra_axonal.weight = model_parameters.extra_axonal.proton_density;
    
elseif (~model_parameters.include_proton_density && model_parameters.include_T1_effect)
    
    E1 = exp(-model_parameters.TR/model_parameters.myelin.T1);
    model_parameters.myelin.weight = sind(model_parameters.flip_angle) * (1 - E1) / (1 - cosd(model_parameters.flip_angle) * E1);
    
    E1 = exp(-model_parameters.TR/model_parameters.intra_axonal.T1);
    model_parameters.intra_axonal.weight = sind(model_parameters.flip_angle) * (1 - E1) / (1 - cosd(model_parameters.flip_angle) * E1);
    
    E1 = exp(-model_parameters.TR/model_parameters.extra_axonal.T1);
    model_parameters.extra_axonal.weight = sind(model_parameters.flip_angle) * (1 - E1) / (1 - cosd(model_parameters.flip_angle) * E1);
    
elseif (model_parameters.include_proton_density && model_parameters.include_T1_effect)
    
    E1 = exp(-model_parameters.TR/model_parameters.myelin.T1);
    model_parameters.myelin.weight = model_parameters.myelin.proton_density * sind(model_parameters.flip_angle) * (1 - E1) / (1 - cosd(model_parameters.flip_angle) * E1);
    
    E1 = exp(-model_parameters.TR/model_parameters.intra_axonal.T1);
    model_parameters.intra_axonal.weight = model_parameters.intra_axonal.proton_density * sind(model_parameters.flip_angle) * (1 - E1) / (1 - cosd(model_parameters.flip_angle) * E1);
    
    E1 = exp(-model_parameters.TR/model_parameters.extra_axonal.T1);
    model_parameters.extra_axonal.weight = model_parameters.extra_axonal.proton_density * sind(model_parameters.flip_angle) * (1 - E1) / (1 - cosd(model_parameters.flip_angle) * E1);
    
end
end