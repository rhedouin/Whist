function [hist, edges] = createHistogramFieldPerturbation(Model, Field, options)
   
    if ~exist('options', 'var')
        options.null = 1;
    end

    if isfield(options, 'new_figure')
        figure;
    else 
        hold on
    end
    
    if isfield(options, 'edges')
        edges.intra_axonal = options.edges.intra_axonal;
        edges.myelin = options.edges.myelin;
        edges.extra_axonal = options.edges.extra_axonal;
    else
        edges.intra_axonal = (-15 : 0.5 : 5);
        edges.myelin = (-10 : 0.5 : 20);
        edges.extra_axonal = (-10 : 0.5 : 10);
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
    hold on
    
    listField = Field(:);
    
    nb_pixels = sum(mask(:));

    [hist.intra_axonal, edges.intra_axonal] = histcounts(listField(Model == 0.5),edges.intra_axonal);
    hist.intra_axonal = hist.intra_axonal/nb_pixels;
    
    [hist.myelin, edges.myelin] = histcounts(listField(Model == 1),edges.myelin);
    hist.myelin = hist.myelin/nb_pixels;
    
    [hist.extra_axonal, edges.extra_axonal] = histcounts(listField(Model == 0),edges.extra_axonal);
    hist.extra_axonal = hist.extra_axonal/nb_pixels;
    
       
    if plot_hist
        plot(edges.intra_axonal(1:end-1),hist.intra_axonal , [line_style 'r'], ...
            'LineWidth',LineWidth)
        
        plot(edges.myelin(1:end-1),hist.myelin , [line_style 'b'], ...
            'LineWidth',LineWidth)
        
        plot(edges.extra_axonal(1:end-1), hist.extra_axonal, [line_style 'g'], ...
            'LineWidth',LineWidth)
    end
    
    leg = legend('intra axonal', 'myelin', 'extra axonal');
    xlabel('Hz')
    title('Frequency histogram')
    
    if isfield(options, 'xlim')
        xlim(options.xlim);
    end
    if isfield(options, 'fontSize')
        set(gca, 'FontSize', options.fontSize);
    end
    if isfield(options, 'fontWeight')
        set(gca, 'FontWeight', options.fontWeight);
    end

end