% Estimate mean parameter values among samples
clear 
close all

parameter_list = {'FVF', 'gRatio', 'T2myel', 'T2out', 'weight', 'xa', 'xi'};
suffix_list = {'_fix_xa_with_normalization', '_fix_xa_without_normalization'}

spinal_cord_list = zeros(2,7);

spinal_cord_list(1,[2 3 5]) = 1;
spinal_cord_list(2,[1 4 5]) = 1;

% parameter_list = {'xi'};

for kpar = 1:length(parameter_list)
    parameter = parameter_list{kpar};
    figure
    it = 0;

    for ksuf = 1:length(suffix_list);

        for kses = 1:2
        
            suffix = suffix_list{ksuf};
            
            it = it+1;
            clear mean_sample mean_all_sample
            
            base_folder  = '/project/3015069.01/derived/Porcine-1/';
            ses = ['ses-mri0' num2str(kses)]
            input_folder = [base_folder ses '/concatenate_signals/'];
            parameter_folder = [input_folder 'parameter_maps/sample_mask' suffix '/'];
            
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
            
            plot(it + ksuf-1, mean_sample_spinal_cord, 'bx')
            hold on
            plot(it + ksuf-1, mean(mean_sample_spinal_cord), 'bx', 'MarkerSize',10, 'LineWidth', 4)
            
            plot(it + ksuf-1, mean_sample_WM, 'gx')
            plot(it + ksuf-1, mean(mean_sample_WM), 'gx', 'MarkerSize',10, 'LineWidth', 4)
            
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
                ylim([0 5])
            elseif  kpar == 6
                ylim([-0.25e-6 0 ])
            elseif  kpar == 7
                ylim([-0.2e-6 0.2e-6 ])
            end
                       
            plot(it + ksuf-1, mean_all_sample, 'rx', 'MarkerSize',10, 'LineWidth', 4)
%             title([suffix ' , ' parameter])
            title(['with norm         ' parameter '       without norm'])

            xlim([0 6])
        end
    end
end