% test read json file
clear
close all
base_folder = '/project/3015069.04/deep_learning/multi_orientations/Porcine2/lowres/';

nb_TE_list = [5, 10, 15, 20, 25, 30];
noise_list = [0, 1, 2];

for k = 1:length(nb_TE_list)
    for l = 1:length(noise_list)
        k
        l
        nb_TE = nb_TE_list(k);
        noise = noise_list(l);
        
        current_folder = [base_folder 'SignalWithNoise' num2str(noise) '_8rep_6orientations_' num2str(nb_TE) 'TE_Porcine2_fix_xa_large_FVF_20_directions_polyfit_cartesian_with_theta_with0reg/'];        
        json_file = [current_folder 'total_history.json'];
        
        fid = fopen(json_file);
        raw = fread(fid,inf);
        str = char(raw');
        current_structure = jsondecode(str);
        current_loss = current_structure.loss;
        current_val_loss = current_structure.val_loss;

        total_last_loss(k, l) = current_loss(end);
        total_last_val_loss(k, l) = current_val_loss(end);
        
        total_min_loss(k, l) = min(current_loss);
        total_min_val_loss(k, l) = min(current_val_loss);
    end
end

colors = linspecer(3)

figure(1)
hold on
for l = 1:length(noise_list)
%     plot(nb_TE_list, total_min_loss(:, l), '--', 'Color', colors(l, :), 'LineWidth', 3);
    h(l) = plot(nb_TE_list, total_min_val_loss(:, l), '-', 'Color', colors(l, :), 'LineWidth', 4);
end

leg = legend(h, '0 %', '1 %', '2 %');
title(leg, 'noise')
ylim([0 0.25])

xlabel('number of TE')
ylabel('val loss')
set(gca, 'FontSize', 32, 'FontWeight','bold')







