function [hist, edges, meanShift, posShift, negShift] = createHistogramFieldPerturbation(Model, Field, options)
   
    if ~exist('options', 'var')
        options.null = 1;
    end

    if options.plot ~= 0
        h = figure('Name', 'Frequency histogram');
        hold on
    end
    
    if isfield(options, 'edges')
        edges = options.edges;
    else
        edges = (-15 : 0.5 : 15);
    end
    
    if isfield(options, 'line_style')
        line_style = options.line_style;
    else
        line_style = '-';
    end
    
    if isfield(options, 'LineWidth')
        LineWidth  = options.LineWidth ;
    else
        LineWidth  = 0.5;
    end
    
    if isfield(options, 'plot')
        plot_hist  = options.plot ;
    else
        plot_hist  = 1;
    end
    
    if isfield(options, 'mask')
        mask  = options.mask ;
        Model(mask == 0) = -1;
    else
        mask = ones(size(Model));
    end
    
    listField = Field(:);
    
    nb_pixels = sum(mask(:));

    [hist.intra_axonal, edges] = histcounts(listField(Model == 0.5),edges);
    hist.intra_axonal = hist.intra_axonal/nb_pixels;
    
    [hist.myelin, edges] = histcounts(listField(Model == 1),edges);
    hist.myelin = hist.myelin/nb_pixels;
    
    [hist.extra_axonal, edges] = histcounts(listField(Model == 0),edges);
    hist.extra_axonal = hist.extra_axonal/nb_pixels;
           
    if plot_hist
        plot(edges(1:end-1),hist.intra_axonal , [line_style 'r'], ...
            'LineWidth',LineWidth)
        
        plot(edges(1:end-1),hist.myelin , [line_style 'b'], ...
            'LineWidth',LineWidth)
        
        plot(edges(1:end-1), hist.extra_axonal, [line_style 'g'], ...
            'LineWidth',LineWidth)
        
        leg = legend('intra axonal', 'myelin', 'extra axonal');
        xlabel('Hz')
        title('Frequency histogram')
    end
    
    
    if isfield(options, 'xlim')
        xlim(options.xlim);
    end
    if isfield(options, 'fontSize')
        set(gca, 'FontSize', options.fontSize);
    end
    if isfield(options, 'fontWeight')
        set(gca, 'FontWeight', options.fontWeight);
    end
    
    meanShift = mean(listField(Model == 1));
    maxShift = max(listField(Model == 1));
    minShift = min(listField(Model == 1));
    
    posShift = maxShift - meanShift;
    negShift = minShift - meanShift;

end