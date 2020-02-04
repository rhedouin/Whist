clear
signalFolder = '/project/3015069.01/model/data/dicos/realisticDico/N400/signals/';
cd(signalFolder) 

temp = h5read('SignalWithNoise1_5rep_fix_rho_without_theta_9orientation.h5py', '/directionsValues');
thetaValues = squeeze(acos(abs(temp(3,:,:,:,:,:,:,:))));

m = matfile('SignalWithNoise1_5rep_fix_rho_without_theta_9orientation_theta_output_parameter_scale.h5py','Writable',true);
m.thetaValues = thetaValues;
keyboard;



directionsValues = permute(dirWithRotation, [5 1 2 3 4]);

directionsValues_bis = repmat(directionsValues,1,1,1,1,10,4,7,5,8);

toto = '/project/3015069.01/model/data/dicos/realisticDico/N400/deepLearning/SignalWithNoise1_5rep_fix_rho_without_theta_9orientation_theta_output/SignalWithNoise1_5rep_fix_rho_without_theta_9orientation_theta_output_parameter_scale.h5py'
h5disp(toto)

s = size(SignalComponent)
for j = 1:10
    plot(angle(SignalComponent(randi(s(1)), randi(s(2)), randi(s(3)), randi(s(4))).Myelin))
    hold on
end


toto = SignalComponent.Myelin;