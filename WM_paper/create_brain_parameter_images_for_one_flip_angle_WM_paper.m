% Create figure from brain for WM paper
clear
close all

concatenate_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/concatenate_signals_9_orientations/';
parameter_folder = [concatenate_folder 'parameter_maps/noise4/'];

parameter_list = {'FVF', 'gRatio', 'T2Myelin', 'T2IntraExtraAxonal', 'weight', 'xiMyelin'};

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_ses-03_gre_orientation-4_brain_mask_all_register_coronal.nii.gz'));

fa = 'fa-35';
flip_angle_folder = [parameter_folder fa '/'];

for k = 1:length(parameter_list)
    parameter = parameter_list{k};
    parameter_map.(parameter) = load_nii_img_only([flip_angle_folder 'BrainSample-2_ses-03_' parameter '_' fa '_20_directions_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    parameter_map.(parameter)(mask==0) = 10000;
    
    if ((k == 3) || (k == 4))
        parameter_map.(parameter) = parameter_map.(parameter)*1e3;
    end
    
    parameter_map_permute.(parameter) = permute(parameter_map.(parameter), [3 2 1]);
    parameter_map_permute.(parameter) = parameter_map_permute.(parameter)(end:-1:1, end:-1:1, :);
end

xslice = 0;
yslice = 128 - 65;
zslice = 0;

edge_x = 0;
edge_y = 0;
edge_z = 20;

for k = 1:length(parameter_list)
    parameter = parameter_list{k};
    
    figure('Name',parameter)
    imagesc(squeeze(parameter_map_permute.(parameter)(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));
    
    colormap('gray');
    cb(k) = colorbar('south');
    x1=get(gca,'position');
    x=get(cb(k),'Position');
    x(1) = 0.2;
    x(2) = 0.01;
    x(3) = 0.6;
    x(4) = 0.03;
    
    set(cb(k),'Position',x)
    set(gca,'position',x1)
    if k == 1
        caxis([0 0.6])
    elseif k == 2
        caxis([0.5 0.85])
    elseif k == 3
        caxis([0 20])
    elseif k == 4
        caxis([0 100])
    elseif k == 5
        caxis([0 3])
    elseif k == 6
        caxis([-0.2 0.2])
    end
    
    axis off
    %     title(parameter)
    set(gca, 'FontSize', 20)
    return;

end

t1 = load_nii_img_only([concatenate_folder 'BrainSample-2_ses-03_mp2rage_orientation-4_t1_brain_2_BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1_brain_without_external_csf_register.nii.gz']);
t1(mask==0) = 10000;

t1_permute = permute(t1, [3 2 1]);
t1_permute = t1_permute(end:-1:1, end:-1:1, :);

    
magn_TE1 = load_nii_img_only([concatenate_folder 'BrainSample-2_ses-03_gre_orientation-4_fa-20_magn_echo-1_without_external_csf_register.nii.gz']);
magn_TE1(mask==0) = 10000;

magn_TE1_permute = permute(magn_TE1, [3 2 1]);
magn_TE1_permute = magn_TE1_permute(end:-1:1, end:-1:1, :);


figure('Name','magn TE1')
imagesc(squeeze(magn_TE1_permute(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));
    
colormap('gray');
% cb = colorbar;
% x1=get(gca,'position');
% x=get(cb,'Position');
% x(2) = 0.2;
% x(4) = 0.6;
% x(3) = 0.03;
% set(cb,'Position',x)
% set(gca,'position',x1)

caxis([800 1800])

axis off
set(gca, 'FontSize', 20)


figure('Name', 't1')
% imagesc(squeeze(t1_permute(1+edge_x:end-edge_x, 1+edge_y:end-edge_y, zslice)));
    imagesc(squeeze(t1_permute(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));

colormap('gray');
cb(k) = colorbar('south');
x1=get(gca,'position');
x=get(cb(k),'Position');
x(1) = 0.2;
x(2) = 0.01;
x(3) = 0.6;
x(4) = 0.03;

set(cb(k),'Position',x)
set(gca,'position',x1)

caxis([0.2 0.7])

axis off
set(gca, 'FontSize', 20)








