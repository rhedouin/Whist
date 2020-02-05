function [axon_collection, FVF_current] = removeAxons(axon_collection, FVF_expected, tol, mask, plot_model)

[~, ~, FVF_current] = createModelFromData(axon_collection, mask, plot_model);
iter_tolerance = 0;
iter_total = 0;

axon_collection_old = axon_collection;

while or((FVF_current < FVF_expected - tol),(FVF_current > FVF_expected))
    iter_total = iter_total + 1;
    display(['model with ' num2str(length(axon_collection)) ' axons, current FVF : ' num2str(FVF_current)]);
    
    if (FVF_current > FVF_expected)
        axon_collection_old = axon_collection;
    elseif (FVF_current < FVF_expected - tol)
        axon_collection = axon_collection_old;

        iter_tolerance = iter_tolerance +1;
        if (iter_tolerance > 50)
            iter_tolerance = 0;
            tol = tol*2;
            display(['double the tolerance : ' num2str(tol);]);
        end 
    else
        error('behavior not expected');
    end
    
    n = randi(length(axon_collection));
    
    axon_collection(n) = [];   
    if (mod(iter_total, 5) == 0)
        [~, ~, FVF_current] = createModelFromData(axon_collection, mask, plot_model);
        pause(0.01);
    else
        [~, ~, FVF_current] = createModelFromData(axon_collection, mask, 0);
    end
    
end
display(['Final model with ' num2str(length(axon_collection)) ' axons, current FVF : ' num2str(FVF_current)]);

end