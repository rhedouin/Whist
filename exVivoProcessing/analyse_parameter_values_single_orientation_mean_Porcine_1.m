% Estimate mean parameter values among samples
clear 
close all

parameter_list = {'FVF', 'gRatio', 'T2myel', 'T2out', 'weight', 'xi', 'dispersion'};
suffix_list = {'_fix_xa'}

spinal_cord_list = zeros(2,7);

spinal_cord_list(1,[2 3 5]) = 1;
spinal_cord_list(2,[1 4 5]) = 1;

for kpar = 1:length(parameter_list)
    parameter = parameter_list{kpar};
    
    figure
    it = 0;
    for kses = 1:2
        
        suffix = '_fix_xa_mean';
        
        it = it+1;
        clear mean_sample mean_all_sample
        
        base_folder  = '/project/3015069.01/derived/Porcine-1/';
        ses = ['ses-mri0' num2str(kses)]
        input_folder = [base_folder ses '/concatenate_signals/'];
        parameter_folder = [input_folder 'parameter_maps/single_orientation/mean/'];
        
        prefix = ['Porcine-1_' ses '_'];
        cd(input_folder)
        
        sample_mask = load_nii_img_only([input_folder 'Porcine-1_' ses '_anat_ref_samples.nii.gz']);
        parameter_map = load_nii_img_only([parameter_folder 'Porcine-1_' ses '_' parameter '_sample_mask' suffix '.nii.gz']);
        
        nb_sample = length(unique(sample_mask)) - 1;
        mean_sample_spinal_cord = [];
        mean_sample_WM = [];
        
        for k = 1:nb_sample
            sample =  single(sample_mask == k);
            sample_nb_pixel = sum(sample, 'all');
            
            if  spinal_cord_list(kses, k) == 1
                mean_sample_spinal_cord(end+1) = sum(parameter_map.*sample, 'all') / sample_nb_pixel;
            else
                mean_sample_WM(end+1) = sum(parameter_map.*sample, 'all') / sample_nb_pixel;
            end
        end
        
        plot(it, mean_sample_spinal_cord, 'bx')
        hold on
        plot(it, mean(mean_sample_spinal_cord), 'bx', 'MarkerSize',10, 'LineWidth', 4)
        
        plot(it, mean_sample_WM, 'gx')
        plot(it, mean(mean_sample_WM), 'gx', 'MarkerSize',10, 'LineWidth', 4)
        
        mean_all_sample = mean([mean_sample_spinal_cord mean_sample_WM]);
        
        if kpar == 1
            ylim([0 1])
        elseif  kpar == 2
            ylim([0 1])
        elseif  kpar == 3
            ylim([0 0.02])
        elseif  kpar == 4
            ylim([0 0.12])
        elseif  kpar == 5
            ylim([0 4])
        elseif  kpar == 6
            ylim([-0.2 0.2 ])
        elseif  kpar == 7
            ylim([0 0.4 ])
        end
        
        plot(it, mean_all_sample, 'rx', 'MarkerSize',10, 'LineWidth', 4)
        title(['single orientation mean: ' parameter])
        xlim([0 3])
    end
end