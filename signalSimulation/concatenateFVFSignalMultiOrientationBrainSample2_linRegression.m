clear
close all

N = 400;
base_folder = '/project/3015069.04/';
signal_folder = [base_folder 'signal_components/multi_orientations/BrainSample2/without_dispersion/'];
dico_folder = [base_folder 'dictionaries/multi_orientations/BrainSample2/'];

num=1;
suffix = ['_train' num2str(num)];
FVF = 40;
FVF_folder = [signal_folder 'FVF' num2str(FVF) '_N' num2str(N) suffix '/'];

load([FVF_folder 'Signal_FVF' num2str(FVF) suffix '.mat']);

FVFRange = [10 20 30 40 50 60 70 80 85];
T2myelRange = [4 8 16 32 64]*1e-3;
T2outRange = [15 30 60 120 340 680]*1e-3;
weightRange = [0.3 0.6 0.9 1.2 1.8 2.4 4.8];

time = linspace(2.15,35.7,12)*1e-3;

lFVF = length(FVFRange);
lgRatio = length(gRatioRange);
lXi = length(xiRange);
lXa = length(xaRange);
lDir = size(fiber_directions,1);
lT2myel = length(T2myelRange);
lT2out = length(T2outRange);
lWeight = length(weightRange); 

nb_orientation = 9;

include_theta = 1;
include_phase = 1;

length_one_orientation = 25;

lVec = nb_orientation * length_one_orientation;

nb_replica = 1;

SignalValues=zeros(lVec, lFVF, lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');

FVFDicoTemp = zeros(lFVF, lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
gRatioDicoTemp = zeros(lFVF, lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
xiDicoTemp = zeros(lFVF, lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
xaDicoTemp = zeros(lFVF, lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
T2myelDicoTemp = zeros(lFVF, lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
T2outDicoTemp = zeros(lFVF, lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
weightDicoTemp = zeros(lFVF, lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
directionsDicoTemp = zeros(3,lFVF, lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');
thetaDicoTemp = zeros(lFVF, lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, nb_replica,'single');

it = 0;
noise = 0;

for nb = 1:nb_replica
    nb
    suffix = ['_train' num2str(nb)];
    
    for k = 1:lFVF
        FVF = FVFRange(k)
        FVF_folder = [signal_folder 'FVF' num2str(FVF) '_N' num2str(N) suffix '/'];
        
        load([FVF_folder 'Signal_FVF' num2str(FVF) suffix '.mat']);
        tic()
        for l = 1:lgRatio
            gRatio = gRatioRange(l);
            for m = 1:lXi
                xi = xiRange(m);
                for n = 1:lXa
                    xa = xaRange(n);
                    for o = 1:lDir
                        for p = 1:lT2myel
                            T2myel = T2myelRange(p);
                            for q = 1:lT2out
                                T2out = T2outRange(q);
                                for r = 1:lWeight
                                    weight = weightRange(r);
                                    tempSignalVector = zeros(lVec, 1);                                    

                                    for rot = 1:nb_orientation
                                        
                                        Signal0 = SignalComponent(l,m,n,o,rot);
                                        
                                        SignalWithR2 = exp(-time/T2out).*Signal0.Axon + weight*exp(-time/T2myel).*Signal0.Myelin + exp(-time/T2out).*Signal0.Extra;
                                        SignalWithR2 = SignalWithR2 / abs(SignalWithR2(1));
                                        SignalWithR2 = SignalWithR2 + noise*(randn(1,12) + 1i*randn(1,12));
                                        SignalWithR2 = SignalWithR2 / abs(SignalWithR2(1));
                                        
                                        
                                        if include_theta
                                            theta = acos(abs(Signal0.current_dir(3)));
                                            theta_noisy = theta + noise*randn(1);
                                        else
                                            theta_noisy = [];
                                        end
                                        
                                        SignalWithR2 = SignalWithR2 / abs(SignalWithR2(1));
                                        
                                        norm_magn =  abs(SignalWithR2);
                                        phase = angle(SignalWithR2);                                                                                                                           
      
                                        norm_phase = phase(3:end) - phase(1) - (phase(2)-phase(1))*(time(3:end) - time(1))/(time(2) - time(1));

                                        poly_coeff = polyfit(time, phase, 1);
                                        norm_phase_poly = phase - (time*poly_coeff(1) + poly_coeff(2));
                                            
                                        norm_complex_poly = norm_magn.*exp(1i*norm_phase_poly);
                                        norm_magn_poly = abs(norm_complex_poly);
                                        norm_phase_poly = angle(norm_complex_poly);
                                        real_magn_poly = real(norm_complex_poly);
                                        imag_magn_poly = imag(norm_complex_poly);

                                        tempSignalDir = [theta_noisy real_magn_poly imag_magn_poly];  
           
                                        tempSignalVector((rot-1)*length_one_orientation + 1 : rot*length_one_orientation) = tempSignalDir;                                             
                                    end
                                                                        
                                    main_direction = SignalComponent(l,m,n,o,4).current_dir;
                                    
                                    SignalValues(:,k,l,m,n,o,p,q,r,nb) = tempSignalVector;
                                                                       
                                    FVFDicoTemp(k,l,m,n,o,p,q,r,nb) = Signal0.FVF;
                                    gRatioDicoTemp(k,l,m,n,o,p,q,r,nb) = Signal0.gRatio;
                                    xiDicoTemp(k,l,m,n,o,p,q,r,nb) = xi;
                                    xaDicoTemp(k,l,m,n,o,p,q,r,nb) = xa;
                                    T2myelDicoTemp(k,l,m,n,o,p,q,r,nb) = T2myel;
                                    T2outDicoTemp(k,l,m,n,o,p,q,r,nb) = T2out;
                                    weightDicoTemp(k,l,m,n,o,p,q,r,nb) = weight;
                                    directionsDicoTemp(:,k,l,m,n,o,p,q,r,nb) = main_direction;
                                    thetaDicoTemp(k,l,m,n,o,p,q,r,nb) = acos(abs(main_direction(3)));
                                end
                            end
                        end
                    end
                end
            end
        end
        toc()
    end
end


SignalValues = squeeze(SignalValues);

FVFValues = squeeze(FVFDicoTemp);
gRatioValues = squeeze(gRatioDicoTemp);
xiValues = squeeze(xiDicoTemp);
xaValues = squeeze(xaDicoTemp);
T2myelValues = squeeze(T2myelDicoTemp);
T2outValues = squeeze(T2outDicoTemp);
weightValues = squeeze(weightDicoTemp);
directionsValues = squeeze(directionsDicoTemp);
thetaValues = squeeze(thetaDicoTemp);

FVF = FVFRange/100;
gRatio = gRatioRange/100;

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

if noise == 0.005
    prefix_name = 'SignalWithNoise05';
else
    prefix_name = ['SignalWithNoise' num2str(100*noise)];
end

base_name = [prefix_name '_' num2str(nb_replica) 'rep_' num2str(nb_orientation) 'orientations_BrainSample2_fix_xa_log_scale_cartesian_polyfit'];
signal_name = [base_name '.h5py'];

save([dico_folder signal_name], 'SignalValues', ...
    'FVFRange', 'gRatioRange', 'xiRange', 'xaRange', 'T2myelRange','T2outRange','weightRange', 'FVFValues', ...
    'gRatioValues', 'xiValues', 'xaValues', 'T2myelValues', 'T2outValues', 'directionsValues', 'thetaValues', ...
    'weightValues', 'nb_replica', 'time', 'infoDico', 'infoSignal', 'sphere_rotations', '-v7.3')

toc()






