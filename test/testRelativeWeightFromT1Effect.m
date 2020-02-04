% Test signal according to different  flip angle 
clear
close all

rho = 0.5;
TR = 55e-3;
T1out = 1.5; 
T1myel = 300 *1e-3;

flip_angle = 90;

ratio = computeRelativeWeightFromT1Effect(T1out, T1myel, rho, TR, flip_angle)
keyboard;

T1outRange = linspace(1, 2.5, 5); 
T1myelRange = linspace(100, 500, 5) *1e-3;

for k = 1:length(T1outRange)
    for l = 1:length(T1myelRange)
        T1out = T1outRange(k);
        T1myel = T1myelRange(l);
    
        ratio(k,l) = computeRelativeWeightFromT1Effect(T1out, T1myel, rho, TR, flip_angle)
    end
end
