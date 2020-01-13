function [axon_collection, FVF_current] = removeAxons(axon_collection, FVF_expected, tol, dims)
plot = 0;
[~, ~, FVF_current] = createModelFromData(axon_collection, dims, plot);
iter_tolerance = 0;
iter_total = 0;

while or((FVF_current < FVF_expected - tol),(FVF_current > FVF_expected))
    iter_total = iter_total + 1;
    display(['model with ' num2str(length(axon_collection)) ' axons, current FVF : ' num2str(FVF_current)]);
    
    if (FVF_current > FVF_expected)
        axon_collection_old = axon_collection;
    elseif (FVF_current < FVF_expected - tol)
        keyboard;
        axon_collection = axon_collection_old;

        iter_tolerance = iter_tolerance +1;
        if (iter_tolerance > 50)
            iter_tolerance = 0;
            tol = tol*2;
        end 
    else
        error('behavior not expected');
    end
    
    n = randi(length(axon_collection));
    
    axon_collection(n) = [];   
    [~, ~, FVF_current] = createModelFromData(axon_collection, dims, plot);

end
display(['Final model with ' num2str(length(axon_collection)) ' axons, current FVF : ' num2str(FVF_current)]);

end