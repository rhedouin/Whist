% Create figure from brain for WM paper
clear
close all

base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/';
ref_folder = [base_folder 'references/'];
concatenate_folder = [base_folder 'concatenate_signals_9_orientations/'];

parameter_list = {'FVF', 'gRatio', 'R2Myelin', 'R2IntraExtraAxonal', 'weight', 'xiMyelin'};

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_mask_sagittal56.nii.gz'));

fa = 'fa-35';
flip_angle_folder = [concatenate_folder 'parameter_maps/BrainSample2LorentzianCorrection/' fa '/'];

for k = 1:length(parameter_list)
    parameter = parameter_list{k};
    parameter_map.(parameter) = load_nii_img_only([flip_angle_folder 'BrainSample-2_ses-03_' parameter '_' fa '_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    
    parameter_map.(parameter)(mask==0) = 10000;

    parameter_map_permute.(parameter) = permute(parameter_map.(parameter), [3 2 1]);
    parameter_map_permute.(parameter) = parameter_map_permute.(parameter)(end:-1:1, end:-1:1, :);
end

xslice = 0;
yslice = 128-72;
zslice = 56;

edge_x = 0;
edge_y = 10;
edge_z = 10;

for k = 1:length(parameter_list)
    parameter = parameter_list{k};
    
    figure('Name',parameter)
    imagesc(squeeze(parameter_map_permute.(parameter)(1+edge_x:end-edge_x, 1+edge_y:end-edge_y, zslice)));
    colormap('gray');

    if k == 1
        caxis([0 0.8])
    elseif k == 2
        caxis([0.5 0.85])
    elseif k == 3
        caxis([30 120])
    elseif k == 4
        caxis([10 60])
    elseif k == 5
        caxis([0.5 2.5])
    elseif k == 6
        caxis([-0.2 0.2])
    end
    
    axis off
    set(gca, 'FontSize', 20)

end





