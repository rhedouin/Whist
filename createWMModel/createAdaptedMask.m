function mask = createAdaptedMask(axonCollection, dims)

    one_mask = ones(dims);
    [~, ~, FVF, ~] = createModelFromData(axonCollection, one_mask, 0);

    safety_coeff = 0.9;
    
    mask_side = (sqrt(FVF * safety_coeff) * dims);
    
    mask = zeros(dims);
    mask(round((dims(1) - mask_side)/2) : round((dims(1) + mask_side)/2), round((dims(2) - mask_side)/2) : round((dims(2) + mask_side)/2)) = 1;
end