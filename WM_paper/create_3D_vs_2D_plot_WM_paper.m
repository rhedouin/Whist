clear
% close all

base_model = '/project/3015069.04/data/2DRM/';
cd(base_model)
load('distance_l2_2D_vs_3D_disp004_new_mask_100_models.mat')

d = sqrt(d);

mean_d = squeeze(mean(d,  3));
min_d = squeeze(min(d,[],  3));
max_d = squeeze(max(d,[],  3));

pos = max_d - mean_d;
neg = min_d - mean_d;

st = 4;
mean_d2 = mean_d(st:end, :);
pos2 = pos(st:end, :);
neg2 = neg(st:end, :);

x = (0:15:90);
x2 = x(st:end);

figure;

for k = 1:size(d,2)
    errorbar(x, mean_d(:,k), neg(:,k), pos(:,k), 'LineWidth', 1);
    hold on
end

% title('l1 histogramm error between 3D and 2D models')
ylabel('m(F_{2D}, F_{3D})')
xlabel('\theta')

leg = legend('0', '0.2', '0.4', '0.6');
leg.NumColumns = 2;

title(leg, 'dispersion')
set(gca, 'FontSize', 14, 'FontWeight','bold' )


axes('position',[.45 .4 .4 .3])

for k = 1:size(d,2)
    errorbar(x2, mean_d2(:,k), neg2(:,k), pos2(:,k), 'LineWidth', 1);
    hold on
end
set(gca, 'FontSize', 10, 'FontWeight','bold' )

return;

