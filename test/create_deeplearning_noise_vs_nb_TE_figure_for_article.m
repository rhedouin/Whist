clear
close all

nb_TE = [5, 10, 15, 20];

val_lost_noise05 = [0.19113561465740203, 0.1549881464799245, 0.13805423504511516, 0.1304800900220871];
 
            
plot(nb_TE, val_lost_noise05, '-', 'LineWidth', 3);
hold on

val_lost_noise1 = [0.2016702551762263, 0.1659479938030243, 0.1487193856239319, 0.14011629763444264];
                      
plot(nb_TE, val_lost_noise1, '-', 'LineWidth', 3);

val_lost_noise2 = [0.2226824606180191, 0.1840070509036382, 0.16397175829410554, 0.1536830891529719];

plot(nb_TE, val_lost_noise2, '-', 'LineWidth', 3);

xlabel('number of TE')
ylabel('val loss')
title('After 10 epochs')
leg = legend('0.5%', '1%', '2%');
title(leg, 'noise')
set(gca, 'FontSize', 24, 'FontWeight','bold')

% ylim([0 0.09])
    