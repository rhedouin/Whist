function [hist edges] = createHistogramFieldPerturbation(Model, Field, options)
   
    if ~exist('options', 'var')
        options.null = 1;
    end

    if isfield(options, 'new_figure')
        figure;
    else 
        hold on
    end
    
    if isfield(options, 'edges')
        edges.intra = options.edges.intra;
        edges.myelin = options.edges.myelin;
        edges.extra = options.edges.extra;
    else
        edges.intra = [-15 : 0.5 : 5];
        edges.myelin = [-10 : 0.5 : 20];
        edges.extra = [-10 : 0.5 : 10];
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
    
    colorValues = linspecer(3);

    
    nb_pixels = sum(mask(:));

    [hist.intra, edges.intra] = histcounts(listField(Model == 0.5),edges.intra);
    hist.intra = hist.intra/nb_pixels;
    
    [hist.myelin, edges.myelin] = histcounts(listField(Model == 1),edges.myelin);
    hist.myelin = hist.myelin/nb_pixels;
    
    [hist.extra, edges.extra] = histcounts(listField(Model == 0),edges.extra);
    hist.extra = hist.extra/nb_pixels;
    
       
    if plot_hist
        plot(edges.intra(1:end-1),hist.intra , line_style, 'color', colorValues(1,:), ...
            'LineWidth',LineWidth)
        
        plot(edges.myelin(1:end-1),hist.myelin , line_style, 'color', colorValues(2,:), ...
            'LineWidth',LineWidth)
        
        plot(edges.extra(1:end-1), hist.extra, line_style, 'color', colorValues(3,:), ...
            'LineWidth',LineWidth)
    end
    
    leg = legend('intra', 'myelin', 'extra')
    xlabel('Hz')
    
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