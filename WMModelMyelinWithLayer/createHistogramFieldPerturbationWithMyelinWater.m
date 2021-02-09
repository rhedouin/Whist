function [hist, meanShift, posShift, negShift] = createHistogramFieldPerturbationWithMyelinWater(Model, Field, options)
   
    if ~exist('options', 'var')
        options.null = 1;
    end

    if (~isfield(options, 'create_figure') || options.create_figure == 1)
        h = figure('Name', 'Frequency histogram');
        hold on
    end
        
    if ~isfield(options, 'edges')
        options.edges = (-15 : 0.5 : 15);
    end
    
    if ~isfield(options, 'line_style')
        options.line_style = '-';
    end
    
    if ~isfield(options, 'LineWidth')
        options.LineWidth = 0.5;
    end
    
    if ~isfield(options, 'plot')
        options.plot = 1;
    end
    
    if  isfield(options, 'mask')
        Model(options.mask == 0) = -1;
    else
        options.mask = ones(size(Model));
    end
    
    listField = Field(:);
    
    nb_pixels = sum(options.mask(:));

    [hist.intra_axonal, ~] = histcounts(listField(Model == 0.5), options.edges);
    hist.intra_axonal = hist.intra_axonal/nb_pixels;
     
%     [hist.myelin_phospholipid, ~] = histcounts(listField(Model == 1), options.edges);
%     hist.myelin_phospholipid = hist.myelin_phospholipid/nb_pixels;  
%     
    [hist.myelin_water, ~] = histcounts(listField(Model == 0.75), options.edges);
    hist.myelin_water = hist.myelin_water/nb_pixels;
    
%     hist.myelin = hist.myelin_water + hist.myelin_phospholipid;
%     
    [hist.extra_axonal, ~] = histcounts(listField(Model == 0), options.edges);
    hist.extra_axonal = hist.extra_axonal/nb_pixels;
           
    if options.plot
        plot(options.edges(1:end-1), hist.extra_axonal, [options.line_style 'g'], ...
            'LineWidth', options.LineWidth)
        plot(options.edges(1:end-1),hist.intra_axonal , [options.line_style 'r'], ...
            'LineWidth', options.LineWidth)
        
%         plot(options.edges(1:end-1),hist.myelin_phospholipid , [options.line_style 'k'], ...
%             'LineWidth', options.LineWidth)
        
        plot(options.edges(1:end-1),hist.myelin_water , [options.line_style 'b'], ...
            'LineWidth', options.LineWidth)
        
%         plot(options.edges(1:end-1),hist.myelin , '--c', ...
%             'LineWidth', options.LineWidth)
        

        
%         leg = legend('intra axonal', 'myelin water', 'extra axonal');
        leg = legend('extra axonal', 'intra axonal', 'myelin water');
        
        xlabel('Hz')
        title('Frequency histogram')
    
    end
    
    if isfield(options, 'xlim')
        xlim(options.xlim);
    end
    if isfield(options, 'ylim')
        ylim(options.ylim);
    end
    
    
    if isfield(options, 'fontSize')
        set(gca, 'FontSize', options.fontSize);
    end
    
    if isfield(options, 'fontWeight')
        set(gca, 'FontWeight', options.fontWeight);
    end

    meanShift = mean(listField(Model == 0.75));
    maxShift = max(listField(Model == 0.75));
    minShift = min(listField(Model == 0.75));
    
    posShift = maxShift - meanShift;
    negShift = minShift - meanShift;
end