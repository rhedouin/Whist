clear
% close all
base_folder = '/project/3015069.04/dictionaries/multi_orientations/Porcine-1/dictionary_part/'

nb_replic = 1;
g_index = 4;
xi_index = 2;
T2_myel_index = 3;
T2_out_index = 3;
weight_index = 3;
nb_dir = 20;

for k = 1:nb_replic
    for dir = 1:nb_dir
        replic_folder = [base_folder 'FVF70_N400_train' num2str(k) '/'];
        signals = h5read([replic_folder 'SignalWithNoise05_FVF70_replic' num2str(k) '_9orientations_Porcine1_fix_xa_polyfit_cartesian_with_theta.h5py'], '/SignalValues');
      
        temp = squeeze(signals(:, g_index, xi_index, dir, T2_myel_index, T2_out_index, weight_index));
        plot(temp)
        hold on
    end
end