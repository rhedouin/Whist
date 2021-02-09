% create variation plot for one model with different noise level WM paper
clear
close all

base_folder = '/project/3015069.04/dictionaries/unique_value/weight2/parameter_maps/';
equal_noise_folder = [base_folder 'adaptive_noise/']
max_noise_folder = [base_folder 'noise4/']

noise_list = {'0', '05', '1', '2', '4'};
parameter_list = {'FVF', 'gRatio', 'xiMyelin', 'T2IntraExtraAxonal', 'T2Myelin', 'weight'};
% parameter_list = {'FVF'};

for k = 1:length(parameter_list)
    parameter = parameter_list{k}
    for l = 1:length(noise_list)
        noise = noise_list{l}
                
        para_map_equal_noise_path = [equal_noise_folder 'cube10_noise' noise '_' parameter '_adaptive_noise.nii.gz'];
        
        para_map_equal_noise = load_nii_img_only(para_map_equal_noise_path);
        para_map_equal_noise = para_map_equal_noise(:);
        
        if (k == 4 ||k == 5)
            para_map_equal_noise = 1000*para_map_equal_noise;
        end
        
        all_para_map_equal_noise.(parameter)(:, l) = para_map_equal_noise;
 
        mean_para_equal_noise(k,l) = mean(para_map_equal_noise);
        std_para_equal_noise(k,l) = std(para_map_equal_noise);

        
        
        para_map_max_noise_path = [max_noise_folder 'cube10_noise' noise '_' parameter '_noise4.nii.gz'];
        
        para_map_max_noise = load_nii_img_only(para_map_max_noise_path);
        para_map_max_noise = para_map_max_noise(:);
        
        if (k == 4 ||k == 5)
            para_map_max_noise = 1000*para_map_max_noise;
        end
        
        all_para_map_max_noise.(parameter)(:, l) = para_map_max_noise;
 
        mean_para_max_noise(k,l) = mean(para_map_max_noise);
        std_para_max_noise(k,l) = std(para_map_max_noise);

    end
end

                                   
dims = size(all_para_map_equal_noise.(parameter));
color = linspecer(dims(2));

original.FVF = 0.7;
original.gRatio = 0.65;
original.T2Myelin = 16;
original.T2IntraExtraAxonal = 60;
original.weight = 2;
original.xiMyelin = -0.1;


for k = 1:length(parameter_list)
    parameter = parameter_list{k}
    y(:, 1, :) = reshape(all_para_map_max_noise.(parameter), [dims(1) 1 dims(2)]);
    y(:, 2, :) = reshape(all_para_map_equal_noise.(parameter), [dims(1) 1 dims(2)]);

    %     y = reshape(all_para_map.(parameter), [dims(1) 1 dims(2)]);
    
    figure;
    h = iosr.statistics.boxPlot({'4% noise','equal noise'},y,...
        'notch',true,...
        'xSeparator',true,...
        'symbolColor','k',...
        'medianColor','k', ...
        'boxcolor',{color(1,:); color(2,:); color(3,:); color(4,:); color(5,:)}, ...
        'groupLabels',{'0 %','0.5 %','1 %','2 %','4 %'}, ...
        'lineWidth', 1.5, ...
        'limit', 'none');
        
    set(gca, 'FontSize',18) 
    
    if k == 1
        ylim([0.3 1])
    elseif k == 2
        ylim([0.5 0.85])
    elseif k == 3
        ylim([-0.3 0.1])
    elseif k == 4
        ylim([20 100])
    elseif k == 5
        ylim([8 24])
    elseif k == 6
        ylim([1 3.5])
    end

    hl(k) = hline(original.(parameter), 'k--')
    hl(k).LineWidth = 1.5;
%     xaxis off
%     title(parameter, 'FontSize', 20)
    
%     ax1 = gca;                   % gca = get current axis
    % ax1.YAxis.Visible = 'off';   % remove y-axis
%     ax1.XAxis.Visible = 'off';   % remove x-axis
end





