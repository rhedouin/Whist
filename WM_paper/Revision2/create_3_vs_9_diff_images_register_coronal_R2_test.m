% Create figure from brain for WM paper
clear
close all

base_folder = '/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/';
orients9_fa35_folder = [base_folder 'concatenate_signals_9_orientations/parameter_maps/BrainSample2LorentzianCorrection/fa-35/'];
orients3_fa35_mean_folder = [base_folder 'concatenate_signals_3_orientations/mean_fa-35/BrainSample2LorentzianCorrection/'];

parameter_list = {'FVF', 'gRatio', 'R2Myelin', 'R2IntraExtraAxonal', 'weight', 'xiMyelin'};
parameter_list = {'FVF'};

mask = single(load_nii_img_only('/project/3015069.01/derived/BrainSample-2/ses-03/gre_renaud/masks/BrainSample-2_mask_coronal_73.nii.gz'));

fa = 'fa-35';

for k = 1:length(parameter_list)
    
    parameter = parameter_list{k};
    parameter_map_orients9.(parameter) = load_nii_img_only([orients9_fa35_folder 'BrainSample-2_ses-03_' parameter '_fa-35_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    parameter_map_orients9.(parameter)(mask==0) = 10000;

    parameter_map_orients9_permute.(parameter) = zeros(82, 128, 128);
    parameter_map_orients9_permute.(parameter)(1:10, :, :) = 10000;
    parameter_map_orients9_permute.(parameter)(11:end, :, :) = permute(parameter_map_orients9.(parameter), [3 2 1]);
    parameter_map_orients9_permute.(parameter) = parameter_map_orients9_permute.(parameter)(end:-1:1, end:-1:1, :);
        
    parameter_map_orients3.(parameter) = load_nii_img_only([orients3_fa35_mean_folder 'BrainSample-2_ses-03_' parameter '_mean_fa-35_polyfit_cartesian_with_theta_noise4_register.nii.gz']);
    parameter_map_orients3.(parameter)(mask==0) = 0;

    parameter_map_orients3_permute.(parameter) = zeros(82, 128, 128);
    parameter_map_orients3_permute.(parameter)(1:10, :, :) = 0;
    parameter_map_orients3_permute.(parameter)(11:end, :, :) = permute(parameter_map_orients3.(parameter), [3 2 1]);
    parameter_map_orients3_permute.(parameter) = parameter_map_orients3_permute.(parameter)(end:-1:1, end:-1:1, :);
    
    diff.(parameter) = abs(parameter_map_orients9_permute.(parameter) - parameter_map_orients3_permute.(parameter));
end

xslice = 0;
yslice = 128-72;
zslice = 0;

edge_x = 0;
edge_y = 0;
edge_z = 10;
 
for k = 1:length(parameter_list)
    parameter = parameter_list{k};
     
    figure('Name',parameter)
    imagesc(squeeze(diff.(parameter)(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));
    colormap('gray');
    cb(k) = colorbar('south');
   
    caxis([0 0.4])
    
    colormap('gray');
    cb(k) = colorbar('south');
    
    parameter = parameter_list{k};
    
    figure('Name',parameter)
    imagesc(squeeze(parameter_map_orients9_permute.(parameter)(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));
    colormap('gray');
    cb(k) = colorbar('south');
    
    caxis([0 0.8])
    
    colormap('gray');
    cb(k) = colorbar('south');
    
    parameter = parameter_list{k};
    
    figure('Name',parameter)
    imagesc(squeeze(parameter_map_orients3_permute.(parameter)(1+edge_x:end-edge_x, yslice, 1+edge_z:end-edge_z)));
    colormap('gray');
    cb(k) = colorbar('south');
    
    caxis([0 0.8])
    
    colormap('gray');
    cb(k) = colorbar('south');

end
return;







