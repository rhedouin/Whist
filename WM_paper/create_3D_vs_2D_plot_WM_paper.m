clear
close all

base_model = '/project/3015069.04/data/2DRM/';
cd(base_model)
load('distance_l2_2D_vs_3D_disp02.mat')

d = sqrt(d);

mean_d = squeeze(mean(d,  3));
min_d = squeeze(min(d,[],  3));
max_d = squeeze(max(d,[],  3));

pos = max_d - mean_d;
neg = min_d - mean_d;

x = [0:15:90];

for k = 1:7
    errorbar(x, mean_d(:,k), neg(:,k), pos(:,k), 'LineWidth', 1);
    hold on
end
title('l1 histogramm error between 3D and 2D models')
ylabel('l1 error')
xlabel('theta')

legend('0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6')
return;

