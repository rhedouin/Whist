%
clear
close all

load('BrainSample2_rotations_ref_2_orientations.mat')
load('20_fiber_orientations_rotations.mat')

C = nchoosek(1:9,3)

for k = 1:size(C,1)
% for k = [66, 47]

    for m = 1:20
        for l=1:3
            display(['k: ' num2str(k) ', l: ' num2str(l) ', m: ' num2str(m)])
            current_rotation = rotations(:,:,C(k,l));
            fiber_rotated = current_rotation*fiber_directions(m,:)';
            theta(l,m) = acos(abs(fiber_rotated(3)));
            all_theta(k,l,m) = theta(l,m);
        end
        theta_sorted(:,m) = sort(theta(:,m));
        diff_theta(m) = theta_sorted(end,m) -  theta_sorted(1,m);
        
    end
%     
%         figure
%         plot(theta_sorted)
        
    avg_diff_theta(k) = mean(diff_theta);
end         

dict_folder = '/project/3015069.04/dictionaries/multi_orientations/BrainSample2/';
dict_path = [dict_folder 'SignalWithNoise4_8rep_9rotations_12TE_BrainSample2_fix_xa_large_FVF_20_directions_polyfit_cartesian_with_theta.h5py'];

xiMyelinValues = h5read(dict_path, '/xiMyelinValues');
xiMyelinRange = h5read(dict_path, '/xiMyelinRange');
xaMyelinValues = h5read(dict_path, '/xaMyelinValues');
xaMyelinRange = h5read(dict_path, '/xaMyelinRange');
weightValues = h5read(dict_path, '/weightValues');
weightRange = h5read(dict_path, '/weightRange');
thetaValues = h5read(dict_path, '/thetaValues');
sphere_rotations = h5read(dict_path, '/sphere_rotations');
nb_replica = h5read(dict_path, '/nb_replica');
infoSignal = h5read(dict_path, '/infoSignal');
infoDico = h5read(dict_path, '/infoDico');
gRatioValues = h5read(dict_path, '/gRatioValues');
gRatioRange = h5read(dict_path, '/gRatioRange');
directionsValues = h5read(dict_path, '/directionsValues');
TE = h5read(dict_path, '/TE');
T2MyelinValues = h5read(dict_path, '/T2MyelinValues');
T2MyelinRange = h5read(dict_path, '/T2MyelinRange');
T2IntraExtraAxonalValues = h5read(dict_path, '/T2IntraExtraAxonalValues');
T2IntraExtraAxonalRange = h5read(dict_path, '/T2IntraExtraAxonalRange');
FVFValues = h5read(dict_path, '/FVFValues');
FVFRange = h5read(dict_path, '/FVFRange');

% SignalValues_save = h5read(dict_path, '/SignalValues');

[a, b] = sort(avg_diff_theta, 'descend');
lVec = 25; 
for k = 1:10
    best_rotation = C(b(k),:)
    list_signal = [];
    for l = 1:3
        list_signal = [list_signal (best_rotation(l) - 1) * lVec + 1 : best_rotation(l) * lVec];
    end
    dict_out_path = [dict_folder 'SignalWithNoise4_8rep_' num2str(best_rotation(1)) num2str(best_rotation(2))  num2str(best_rotation(3))   'rotations_12TE_BrainSample2_fix_xa_large_FVF_20_directions_polyfit_cartesian_with_theta.h5py'];

%     SignalValues = SignalValues_save(list_signal, :, :, :, :, :, :, :, :, :);
%     save(dict_out_path, 'SignalValues', 'xiMyelinValues', 'xiMyelinRange', ...
%         'xaMyelinValues', 'xaMyelinRange', 'weightValues', 'weightRange', ...
%         'thetaValues', 'sphere_rotations', 'nb_replica', 'infoSignal', ...
%         'infoDico', 'gRatioValues', 'gRatioRange', 'directionsValues', ...
%         'TE', 'T2MyelinValues', 'T2MyelinRange', 'T2IntraExtraAxonalValues', ...
%         'T2IntraExtraAxonalRange', 'FVFValues', 'FVFRange', '-v7.3')
%     
    
end


        


