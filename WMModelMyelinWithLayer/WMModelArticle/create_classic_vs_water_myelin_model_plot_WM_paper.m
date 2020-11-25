clear
close all

base_model = '/project/3015069.04/code/Whist/WMModelMyelinWithLayer/WMModelArticle/';
cd(base_model)
% keyboard;
load('simi_signal_N1000.mat')
x = round(linspace(0, 90, size(simi_signal,2)));
for k = 1:size(simi_signal, 1)
    figure;
    for l = 2:size(simi_signal, 3)
        plot(x, simi_signal(k,:,l), 'LineWidth', 2);
        hold on
    end
    ylim([0 0.025])
    xlim([0 100])

    title(['Myelin water proportion = 0.' num2str(k+1)])
    ylabel('RMSE(M1, M2)')
    xlabel('\theta')
    
    leg =  legend('5', '10', '15', '20');
    leg.NumColumns = 2;

    title(leg, {'number of myelin', 'water layers'})
    set(gca, 'FontSize', 14, 'FontWeight','bold', 'box', 'off')
end


    
    