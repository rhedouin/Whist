% test read json file
clear
close all
base_folder = '/project/3015069.04/deep_learning/multi_orientations/Porcine2/lowres/';

noise_list = {'0', '05',  '1', '2', '4'};
nb_TE = 18;

for l = 1:length(noise_list)

    noise = noise_list{l};
    
    current_folder = [base_folder 'SignalWithNoise' noise '_8rep_6orientations_' num2str(nb_TE) 'TE_Porcine2_fix_xa_large_FVF_20_directions_polyfit_cartesian_with_theta_with0reg/'];
    json_file = [current_folder 'total_history.json'];
    
    fid = fopen(json_file);
    raw = fread(fid,inf);
    str = char(raw');
    current_structure = jsondecode(str);
    current_loss = current_structure.loss;
    current_val_loss = current_structure.val_loss;
    
    total_loss(:, l) = current_loss;
    total_val_loss(:, l) = current_val_loss;
end

colors = linspecer(length(noise_list))

figure(1)
hold on
for l = 1:length(noise_list)
    plot(total_loss(:, l), '--', 'Color', colors(l, :), 'LineWidth', 4);
    h(l) = plot(total_val_loss(:, l), '-', 'Color', colors(l, :), 'LineWidth', 4);
end

leg = legend(h, '0 %', '0.5 %','1 %', '2 %', '4 %');
leg.NumColumns = 2;
title(leg, 'noise')
ylim([0 0.25])

xlabel('number of epochs')
ylabel('loss / val loss')
set(gca, 'FontSize', 32, 'FontWeight','bold')







