%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example shows how to simulate a field perturbation from a WM model
% (here a single axon) and the corresponding multi GRE signals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

%%%%%%%%%%%% Load a WM model with a single axon
model_path = '/project/3015069.04/WM_Models/oneAxon/oneAxon.mat';

load(model_path)
dims = size(mask);

% plot the WM model, the fiber volume fraction (FVF) is computed within the
% mask represented by the red rectangle
plot_model = 1;
createModelFromData(oneAxon, mask, plot_model);

%%%%%%%%%%% Set parameters
% Field strength (in Tesla)
B0 = 3;

% Myelin susceptibility (in ppm) relative to intra/extra axonal compartments
% Anisotropic
xa = -0.1;
% Isotropic
xi = -0.1;

% T2 values (in second)
T2_intra_extra = 50*1e-3;
T2_myelin = 15*1e-3;

% TE (in second)
TE = (2:3:59)*1e-3;

% Relative water signal received from myelin compare to intra/extra compartment (see
% computeRelativeWeightFromT1Effect.m)
weight = 1;

% Magnetic field orientation
theta = pi/2;
phi = 0;
field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];

%%%%%%%%%%%%% Simulate the field perturbation from the WM model and the multi GRE signals
[signal_original, field] = simulateSignalFromModel(oneAxon, mask, xa, xi,  T2_intra_extra, T2_myelin, weight, TE, B0, field_direction);

magn_signal = abs(signal_original);
phase_signal = phase(signal_original);

%%%%%%%%%%%%% Plot field and signals
h = figure('Name', 'Field Perturbation');
position = [10 10 990 790];
h.Position = position;    

subplot(221)
imagesc(field)

cb = colorbar;

title(cb,'Hertz');

caxis([-10 10])
title('Field perturbation', 'FontWeight', 'bold')
set(gca, 'FontSize', 15);


subplot(223)
h3 = subplot(223);
axis([-2 2 -2 2 -2 2]);
CameraPosition = [5 0 0];
view([0 90]);

hold on
arrow3([0 0 0], [0 0 1],'--r',0.5,[],0)
arrow3([0 0 0], [0 1 0],'--y',0.5,[],0)
arrow3([0 0 0], [1 0 0],'--b',0.5,[],0)

arrow3([0 0 0],field_direction,'2.5s',1,[],0)

hold off, axis off, camlight left
set(gca,'CameraViewAngle',4)
text(1,0,0,'X'), text(0,1,0,'Y')
text(0,0,1,'Z','VerticalAlignment','bottom',...
    'HorizontalAlignment','center')

text(-0.5,1.5,0,'B0 orientation', 'FontWeight', 'bold', 'FontSize', 12)
set(gca, 'FontSize', 12);


subplot(222)
plot(TE, magn_signal, 'LineWidth', 3)
xlabel('echo time')
ylabel('|S(t)|')
title('Signal magnitude', 'FontWeight', 'bold')
set(gca, 'FontSize', 12);

subplot(224)

plot(TE, phase_signal, 'LineWidth', 3)
xlabel('echo time')
ylabel('phase(S(t))')

title('Signal phase', 'FontWeight', 'bold')
set(gca, 'FontSize', 12);









