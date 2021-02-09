clear
close all

base_model = '/project/3015069.04/data/2DRM/';
cd(base_model)
% keyboard;
load('distance_l2_2D_vs_3D_disp004_new_mask_100_models.mat')
load('signal_and_histo_distance_2D_vs_3D_disp004_new_mask_100_models_paper_parameter_values.mat')

d = diff_histo;
d = simi_signal;
% d = sqrt(d);

mean_d = squeeze(mean(d,  3));
min_d = squeeze(min(d,[],  3));
max_d = squeeze(max(d,[],  3));

pos = max_d - mean_d;
neg = min_d - mean_d;

x = (0:15:90);
x2 = x(st:end);

figure;
for k = 1:size(d,2)
    errorbar(x, mean_d(:,k), neg(:,k), pos(:,k), 'LineWidth', 1);
    hold on
end
ylim([0 0.07])

% title('l1 histogramm error between 3D and 2D models')
ylabel('RMSE(S_{2D}, S_{3D})')
xlabel('\theta')

leg = legend('0', '0.2', '0.4', '0.6');
leg.NumColumns = 2;

title(leg, '2D model dispersion')
set(gca, 'FontSize', 14, 'FontWeight','bold' )


% axes('position',[.45 .4 .4 .3])
% 
% for k = 1:size(d,2)
%     errorbar(x2, mean_d2(:,k), neg2(:,k), pos2(:,k), 'LineWidth', 1);
%     hold on
% end
% set(gca, 'FontSize', 10, 'FontWeight','bold' )

return;

