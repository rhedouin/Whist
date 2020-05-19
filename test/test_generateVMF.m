clear
close all

% kappa_list = [10000 18 9 5.5 3.5 2 1];
% dispersion_list = [0.001 0.1 0.2 0.3 0.4 0.5 0.6];

kappa_list = [5.5 3.5 2 1];
dispersion_list = [0.3 0.4 0.5 0.6];

avrg_orient = [0; 0; 1];

nb_orientations = 100;

for k = 1:length(kappa_list)
    kappa = kappa_list(k);
    dispersion = dispersion_list(k);
    
    generateVMFsample_for_dispersion(avrg_orient, kappa, dispersion, nb_orientations)
%     keyboard;
end