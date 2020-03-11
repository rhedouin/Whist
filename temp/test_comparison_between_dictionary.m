% temp compare signals
clear
% close all

long_TE_folder = '/project/3015069.04/dictionaries/multi_orientations/Porcine2/lowres/fix_xa_large_FVF_20_directions_30_TE/';
normal_TE_folder = '/project/3015069.04/dictionaries/multi_orientations/Porcine2/lowres/fix_xa_large_FVF_20_directions/';

FVF = 80;
num = 3;
long = h5read([long_TE_folder 'SignalWithNoise0_8rep_6orientations_20TE_Porcine2_fix_xa_large_FVF_20_directions_polyfit_cartesian_with_theta.h5py'], '/SignalValues');
normal = h5read([normal_TE_folder 'SignalWithNoise0_8rep_6orientations_18TE_Porcine2_fix_xa_large_FVF_20_directions_polyfit_cartesian_with_theta.h5py'], '/SignalValues');

keyboard;
k = 6;
l = 4;
m = 1;
n = 1;
o = 6;
p = 3;
q = 3;
r = 3;
s = 3;

for s = 1:8
long_concatenate_signal = long(:, k, l, m, n, o, p, q, r, s);
normal_concatenate_signal = normal(:, k, l, m, n, o, p, q, r, s);
figure(2)
subplot(211)

plot(long_concatenate_signal)
hold on 
subplot(212)
hold on
plot(normal_concatenate_signal)

end