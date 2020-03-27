% test read json file
clear
close all
base_folder = '/project/3015069.04/deep_learning/multi_orientations/theorically_good_16_rotations/';

nb_rotation_max = 16;
noise_list = [0, 1, 2];

for k = 1:nb_rotation_max
    for l = 1:length(noise_list)
        k
        l
        noise = noise_list(l);
        
        current_folder = [base_folder 'SignalWithNoise' num2str(noise) '_8rep_12TE_' num2str(k) 'rotations_fix_xa_large_FVF_20_directions_polyfit_cartesian_with_theta/'];        
        json_file = [current_folder 'total_history.json'];
        
        fid = fopen(json_file);
        raw = fread(fid,inf);
        str = char(raw');
        current_structure = jsondecode(str);
        current_loss = current_structure.loss(1:20);
        current_val_loss = current_structure.val_loss(1:20);

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
    plot(total_min_loss(:, l), '--', 'Color', colors(l, :), 'LineWidth', 4);
%     plot(total_last_loss(:, l), '+', 'Color', colors(l, :), 'LineWidth', 3);

    h(l) = plot(total_last_val_loss(:, l), '-', 'Color', colors(l, :), 'LineWidth', 4);
%     h(l) = plot(total_min_val_loss(:, l), '-', 'Color', colors(l, :), 'LineWidth', 4);
end

leg = legend(h, '0 %', '1 %', '2 %');
title(leg, 'noise')
ylim([0 0.25])

xlabel('number of rotations')
ylabel('loss / val loss')
set(gca, 'FontSize', 32, 'FontWeight','bold')







