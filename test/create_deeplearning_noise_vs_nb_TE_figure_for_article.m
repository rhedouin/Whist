clear
close all



val_lost_noise05 = [0.04571451278924942, 0.03944923912286758, 0.03859057043592135, 0.03839211661616961];
 
            
plot(val_lost_noise05, '-', 'LineWidth', 3);
hold on

val_lost_noise_1 = [0.058126393554608025, 0.04446948475639025, 0.04456581537723541, 0.04372276110847791];

                      
plot(val_lost_noise_1, '-', 'LineWidth', 3);



val_lost_noise2 = [0.0818171665430069, 0.05567376306851705, 0.05175923918088277, 0.0524944716334343];
           

plot(val_lost_noise2, '-', 'LineWidth', 3);

xlabel('number of TE')
ylabel('val loss')
title('After 10 epochs')
leg = legend('0.5%', '1%', '2%');
title(leg, 'noise')
set(gca, 'FontSize', 15, 'FontWeight','bold')

ylim([0 0.09])
    