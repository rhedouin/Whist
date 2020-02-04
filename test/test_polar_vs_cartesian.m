signal_path{1} = 'SignalWithNoise05_FVF60_replic1_9orientations_Porcine1_fix_xa_classic_polar_with_theta.h5py';
signal_path{2} = 'SignalWithNoise05_FVF60_replic1_9orientations_Porcine1_fix_xa_polyfit_cartesian_with_theta.h5py';
signal_path{3} = 'SignalWithNoise05_FVF60_replic1_9orientations_Porcine1_fix_xa_polyfit_polar_with_theta.h5py';

SignalValues{1} = h5read(signal_path{1}, '/SignalValues');
SignalValues{2} = h5read(signal_path{2}, '/SignalValues');
SignalValues{3} = h5read(signal_path{3}, '/SignalValues');

a = 7;
b = 4;
c = 1;
d = 4;
e = 3;
f = 3;
% keyboard;
figure

subplot(211)
hold on
plot(SignalValues{1}(:, a, b, c, d, e, f));

subplot(212)
hold on
plot(SignalValues{2}(:, a, b, c, d, e, f));
plot(SignalValues{3}(:, a, b, c, d, e, f));

