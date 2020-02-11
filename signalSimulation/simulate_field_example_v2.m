% example simulate field perturbation and corresponding multi GRE signal
clear
close all

%%%%%%%%%%%% Load a WM model with a single axon
model_path = '/project/3015069.04/WM_Models/oneAxon/oneAxon.mat';

load(model_path)
dims = size(mask);

% plot the WM model, the fiber volume fraction (FVF) is computed within the
% mask represented by the red rectangle
plot_model = 1;
mask = zeros(dims);
mask(30:70, 30:70) = 1;
model = createModelFromData(oneAxon, mask, plot_model);

%%%%%%%%%%% Set parameters
% field strength in Tesla
model_parameters.B0 = 3;

% myelin 
model_parameters.myelin.T2 = 15*1e-3;
model_parameters.myelin.T1 = 500*1e-3;
model_parameters.myelin.proton_density= 0.5; 

model_parameters.myelin.xi = -0.1;  % myelin anisotropic susceptibility (ppm)
model_parameters.myelin.xa = -0.1;  % myelin isotropic susceptibility (ppm)

% Relative water signal received from myelin compare to intra/extra compartment (see
% computeRelativeWeightFromT1Effect.m)

model_parameters.flip_angle = 20;
model_parameters.TR = 60*1e-3;

% intra axonal
model_parameters.intra_axonal.T2 = 50*1e-3;
model_parameters.intra_axonal.T1 = 1.5;
model_parameters.intra_axonal.proton_density= 0.5; 

% extra axonal
model_parameters.extra_axonal.T2 = 50*1e-3;
model_parameters.extra_axonal.T1 = 1.5;
model_parameters.extra_axonal.proton_density= 0.5; 

% TE 
model_parameters.TE = (2:3:59)*1e-3;
model_parameters = computeCompartmentSignalWeight(model_parameters);

% magnetic field orientation
theta = pi/2;
phi = pi/2;
field_direction = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
model_parameters.field_direction = field_direction;

% mask
model_parameters.mask = mask;

%%%%%%%%%%%%% Simulate the field perturbation from the WM model and the multi GRE signals
[signal_original, field] = simulateSignalFromModel_v2(oneAxon, model_parameters);

options.new_figure = 1;
options.mask = mask;

createHistogramFieldPerturbation(model, field, options);

magn_signal = abs(signal_original.total);
phase_signal = phase(signal_original.total);

%%%%%%%%%%%%% Plot field and signals
figure; 

subplot(221)
imagesc(field)

cb = colorbar;

title(cb,'Hertz');

caxis([-10 10])
title('Field perturbation', 'FontWeight', 'bold')
set(gca, 'FontSize', 15);


subplot(223)
h3 = subplot(223)   ;
axis([-2 2 -2 2 -2 2])
CameraPosition = [5 0 0]
view([0 90])

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
plot(magn_signal, 'LineWidth', 3)
xlabel('echo time')
ylabel('|S(t)|')
title('Signal magnitude', 'FontWeight', 'bold')
set(gca, 'FontSize', 12);

subplot(224)

plot(phase_signal, 'LineWidth', 3)
xlabel('echo time')
ylabel('phase(S(t))')

title('Signal phase', 'FontWeight', 'bold')
set(gca, 'FontSize', 12);









