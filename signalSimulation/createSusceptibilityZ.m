function susceptibility_Z = createSusceptibilityZ(total_X, model, model_parameters)

susceptibility_Z = zeros(model_parameters.dims);

ind_myelin = find(model == 1);

[sub_myelin_x, sub_myelin_y] = ind2sub(model_parameters.dims, ind_myelin);

for k = 1:size(sub_myelin_x,1)
    susceptibility_Z(sub_myelin_x(k), sub_myelin_y(k)) = model_parameters.field_direction * squeeze(total_X(sub_myelin_x(k), sub_myelin_y(k),:,:)) * model_parameters.field_direction'; 
end
end