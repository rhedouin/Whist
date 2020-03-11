% temp compare signals
clear
% close all

long_TE_folder = '/project/3015069.04/signal_components/multi_orientations/Porcine2/lowres/fix_xa_large_FVF_20_directions_30_TE/';
normal_TE_folder = '/project/3015069.04/signal_components/multi_orientations/Porcine2/lowres/fix_xa_large_FVF_20_directions/';

FVF = 80;
num = 3;
long = load([long_TE_folder 'FVF' num2str(FVF) '_N400_train' num2str(num) '/' 'Signal_FVF' num2str(FVF) '_train' num2str(num) '.mat']);
normal = load([normal_TE_folder 'FVF' num2str(FVF) '_N400_train' num2str(num) '/' 'Signal_FVF' num2str(FVF) '_train' num2str(num) '.mat']);

k = 5;
l = 3;
m = 1;
n = 1;
o = 6;

long_intra_axonal = long.signal_components(k, l, m, n, o).intra_axonal.signal;
normal_intra_axonal = normal.signal_components(k, l, m, n, o).intra_axonal.signal;
figure(1)
subplot(211)
plot(abs(normal_intra_axonal))
hold on 
plot(abs(long_intra_axonal))

subplot(212)
plot(angle(normal_intra_axonal))
hold on 
plot(angle(long_intra_axonal))

figure(2)
long_extra_axonal = long.signal_components(k, l, m, n, o).extra_axonal.signal;
normal_extra_axonal = normal.signal_components(k, l, m, n, o).extra_axonal.signal;
subplot(211)
plot(abs(normal_extra_axonal))
hold on 
plot(abs(long_extra_axonal))

subplot(212)
plot(angle(normal_extra_axonal))
hold on 
plot(angle(long_extra_axonal))


figure(3)
long_myelin = long.signal_components(k, l, m, n, o).myelin.signal;
normal_myelin = normal.signal_components(k, l, m, n, o).myelin.signal;
subplot(211)
plot(abs(long_myelin))
hold on 
plot(abs(normal_myelin))

subplot(212)
plot(angle(long_myelin))
hold on 
plot(angle(normal_myelin))
