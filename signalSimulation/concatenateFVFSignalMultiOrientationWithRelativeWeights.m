clear

N = 400;
base_folder = '/project/3015069.01/model/data/dicos/';
signal_folder = [base_folder 'signal_components/multi_orientations/Porcine-1/without_dispersion/'];
dico_folder = [base_folder 'dictionaries/multi_orientations/Porcine-1/'];

num=1;
suffix = ['_train' num2str(num)];
FVF = 40;
FVF_folder = [signal_folder 'FVF' num2str(FVF) '_N' num2str(N) suffix '/'];

load([FVF_folder 'Signal_FVF' num2str(FVF) suffix '.mat']);

FVFRange = (40 : 5 : 85)';
T2myelRange = (4 : 4 : 20)*1e-3;
% T2myelRange = 12*1e-3;

T2outRange = (20 : 15 : 80)*1e-3;
weightRange = (0.2 : 0.2 : 2);

lFVF = length(FVFRange);
lgRatio = length(gRatioRange);
lXi = length(xiRange);
lDir = size(SignalComponent,3);
lT2myel = length(T2myelRange);
lT2out = length(T2outRange);
lWeight = length(weightRange); 
thetaRange = acos(squeeze(dirWithRotationValues(1,1,:,4,3)));

nb_orientation = 9;

include_theta = 1;
include_phase = 1;

if (include_theta && include_phase)
    length_one_orientation = 23;
elseif (~include_theta && include_phase)
   length_one_orientation = 22;
elseif (~include_theta && ~include_phase)
   length_one_orientation = 12;
else 
    error('unknown situation')
end

lVec = nb_orientation * length_one_orientation;

time = linspace(2.15,35.7,12)*1e-3;
tic()
toc()

nb_replica = 5;

SignalValues=zeros(lVec, lFVF, lgRatio, lXi, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
FVFDicoTemp = zeros(lFVF, lgRatio, lXi, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
gRatioDicoTemp = zeros(lFVF, lgRatio, lXi, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
xiDicoTemp = zeros(lFVF, lgRatio, lXi, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
T2myelDicoTemp = zeros(lFVF, lgRatio, lXi, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
T2outDicoTemp = zeros(lFVF, lgRatio, lXi, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
weightDicoTemp = zeros(lFVF, lgRatio, lXi, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
directionsDicoTemp = zeros(3,lFVF, lgRatio, lXi, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
thetaDicoTemp = zeros(lFVF, lgRatio, lXi, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');

noise = 0.01;
for nb = 1:nb_replica
    nb
    tic()
    suffix = ['_train' num2str(nb)];
    
    for j = 1:lFVF
        FVF = FVFRange(j)
        FVF_folder = [signal_folder 'FVF' num2str(FVF) '_N' num2str(N) suffix '/'];
        
        load([FVF_folder 'Signal_FVF' num2str(FVF) suffix '.mat']);
        
        for k = 1:lgRatio
            gRatio = gRatioRange(k);
            for l = 1:lXi
                xi = xiRange(l);
                for m = 1:lDir
                    for n = 1:lT2myel
                        T2myel = T2myelRange(n);
                        for o = 1:lT2out
                            T2out = T2outRange(o);
                            for p = 1:lWeight
                                weight = weightRange(p);
                                tempSignalVector = zeros(lVec, 1);
                                
                                for rot = 1:nb_orientation
                                    
                                    Signal0 = SignalComponent(k,l,m,rot);
                                                                        
                                    SignalWithR2 = weight*exp(-time/T2out).*Signal0.Axon + exp(-time/T2myel).*Signal0.Myelin + weight*exp(-time/T2out).*Signal0.Extra;
                                    SignalWithR2 = SignalWithR2 / abs(SignalWithR2(1));
                                    SignalWithR2 = SignalWithR2 + noise*(randn(1,12) + i*randn(1,12));
                                    SignalWithR2 = SignalWithR2 / abs(SignalWithR2(1));
                                    
                                    if include_theta
                                        theta = acos(abs(dirWithRotationValues(k,l,m,rot,3)));
                                        theta_noisy = theta + noise*randn(1);
                                    else
                                        theta_noisy = [];
                                    end
                                    
                                    norm_magn =  abs(SignalWithR2);
                                    
                                    if include_phase
                                        phase = angle(SignalWithR2);
                                        norm_phase = phase(3:end) - phase(1) - (phase(2)-phase(1))*(time(3:end) - time(1))/(time(2) - time(1));
                                    else
                                        norm_phase = [];
                                    end
                                    
                                    tempSignalDir = [theta_noisy norm_magn norm_phase];
                                    tempSignalVector((rot-1)*length_one_orientation + 1 : rot*length_one_orientation) = tempSignalDir;
                                    
                                end
                                
                                main_direction = Signal0(k,l,m,9).current_dir;
                                
                                SignalValues(:,j,k,l,m,n,o,p,nb) = tempSignalVector;
                                
                                FVFDicoTemp(j,k,l,m,n,o,p,nb) = Signal0.FVF;
                                gRatioDicoTemp(j,k,l,m,n,o,p,nb) = Signal0.gRatio;
                                xiDicoTemp(j,k,l,m,n,o,p,nb) = xi;
                                T2myelDicoTemp(j,k,l,m,n,o,p,nb) = T2myel;
                                T2outDicoTemp(j,k,l,m,n,o,p,nb) = T2out;
                                weightDicoTemp(j,k,l,m,n,o,p,nb) = weight;
                                directionsDicoTemp(:,j,k,l,m,n,o,p,nb) = main_direction;
                                thetaDicoTemp(j,k,l,m,n,o,p,nb) = acos(main_direction(3));
                            end
                        end
                    end
                end
            end
        end
    end
    toc()
end
SignalValues = squeeze(SignalValues);

FVFValues = squeeze(FVFDicoTemp);
gRatioValues = squeeze(gRatioDicoTemp);
xiValues = squeeze(xiDicoTemp);
T2myelValues = squeeze(T2myelDicoTemp);
T2outValues = squeeze(T2outDicoTemp);
weightValues = squeeze(weightDicoTemp);
directionsValues = squeeze(directionsDicoTemp);
thetaValues = squeeze(thetaDicoTemp);

infoDico = 'In order, FVF, gRatio, xi, dir, T2myel, T2out, weight';
infoSignal = ['concatenation of ' num2str(nb_orientation) ' rotations each composed of the theta angle (B0 angle) 12 TE signal magnitude and 10 TE normalized angles'];

tic()
if (include_theta && include_phase)
    vector_type_name = 'with_theta';
elseif (~include_theta && include_phase)
    vector_type_name = 'without_theta';
elseif (~include_theta && ~include_phase)
    vector_type_name = 'magn_only';
else 
    error('unknown situation')
end

base_name = ['SignalWithNoise' num2str(100*noise) '_' num2str(nb_replica) 'rep_relative_weight_' vector_type_name '_' num2str(nb_orientation) 'orientations_fix_xa'];
signal_name = [base_name '.h5py'];

save([dico_folder signal_name], 'SignalValues', ...
    'FVFRange', 'gRatioRange', 'xiRange','T2myelRange','T2outRange','weightRange', 'FVFValues', ...
    'gRatioValues', 'xiValues', 'T2myelValues', 'T2outValues', 'directionsValues', 'thetaValues', ...
    'weightValues','time', 'infoDico', 'infoSignal', 'sphere_rotations', '-v7.3')

FVF = FVFRange/100;
gRatio = gRatioRange/100;

FVF = max(abs(FVF));
gRatio = max(abs(gRatio));
xi = max(abs(xiRange));
theta = max(thetaRange(:));
T2out = max(abs(T2outRange));
T2myel = max(abs(T2myelRange));
weight = max(abs(weightRange));

parameter_name = [base_name '_parameter_scale.h5py'];

save([dico_folder parameter_name], 'weight', 'FVF', 'gRatio', 'xi', 'T2out', 'T2myel', 'theta', 'nb_replica', '-v7.3')

toc()








