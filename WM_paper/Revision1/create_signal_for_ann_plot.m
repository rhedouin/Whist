close all

figure; 
x = 176 : 225;
% plot(signals(1:55,6,5,2,1,5,4,4,5,3), 'r', 'LineWidth', 2)
plot(x, signals(176:225,6,5,2,1,5,4,4,5,3), 'r', 'LineWidth', 3)
xticks([210 225])
hold on
vline(201, '--k')
% vline(50, '--k')

ylim([-0.25 1.6])
set(gca, 'FontSize', 20, 'FontWeight','bold', 'box', 'off')

figure; 
x = 1 : 55;
plot(x, signals(1:55,6,5,2,1,5,4,4,5,3), 'r', 'LineWidth', 3)
% xticks([210 225])
hold on
vline(25, '--k')
vline(50, '--k')

ylim([-0.25 1.6])
set(gca, 'FontSize', 20, 'FontWeight','bold', 'box', 'off')

