function out = createDictionaryPartSignalMultiModalities(signal_path, output_folder, T2myelRange, T2outRange, weightRange, time, noise, FVF, nb_orientations, replic, dict_options)
load(signal_path);

lgRatio = length(gRatioRange);
lXi = length(xiRange);
lXa = length(xaRange);
lDir = size(fiber_directions,1);
lT2myel = length(T2myelRange);
lT2out = length(T2outRange);
lWeight = length(weightRange); 

nb_TE = length(time);

length_one_orientation = 2*nb_TE + 1;
lVec = nb_orientations * length_one_orientation;

coordinate_total_list = fieldnames(dict_options.coordinate);
it = 0;

for kCoordinate = 1:length(coordinate_total_list)
    if dict_options.coordinate.(coordinate_total_list{kCoordinate}) == 1
        it = it + 1;
        coordinate_list{it} = coordinate_total_list{kCoordinate};
        SignalValuesMultiCoordinates.(coordinate_list{it}) = zeros(lVec, lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, 'single');
    end
end
lCoordinate = length(coordinate_list);

gRatioDicoTemp = zeros(lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, 'single');
xiDicoTemp = zeros(lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, 'single');
xaDicoTemp = zeros(lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, 'single');
T2myelDicoTemp = zeros(lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, 'single');
T2outDicoTemp = zeros(lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, 'single');
weightDicoTemp = zeros(lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, 'single');
directionsDicoTemp = zeros(3,lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, 'single');
thetaDicoTemp = zeros(lgRatio, lXi, lXa, lDir, lT2myel, lT2out, lWeight, 'single');

tic()
for l = 1:lgRatio
    gRatio = gRatioRange(l)
    for m = 1:lXi
        xi = xiRange(m)
        for n = 1:lXa
            xa = xaRange(n)
            tic()
            for o = 1:lDir
                for p = 1:lT2myel
                    T2myel = T2myelRange(p);
                    for q = 1:lT2out
                        T2out = T2outRange(q);
                        for r = 1:lWeight
                            weight = weightRange(r);
                            
                            for kCoordinate = 1:lCoordinate
                                coordinate = coordinate_list{kCoordinate};
                                tempSignalVector.(coordinate) = zeros(lVec, 1);
                            end
                                                        
                            for rot = 1:nb_orientations
                                
                                Signal0 = SignalComponent(l,m,n,o,rot);

                                SignalWithR2 = exp(-time/T2out).*Signal0.Axon + weight*exp(-time/T2myel).*Signal0.Myelin + exp(-time/T2out).*Signal0.Extra;
                                SignalWithR2 = SignalWithR2 / abs(SignalWithR2(1));
                                SignalWithR2 = SignalWithR2 + noise*(randn(1,nb_TE) + 1i*randn(1,nb_TE));
                                SignalWithR2 = SignalWithR2 / abs(SignalWithR2(1));
                                
                                if dict_options.include_theta
                                    theta = acos(abs(Signal0.current_dir(3)));
                                    theta_noisy = theta + noise*randn(1);
                                else
                                    theta_noisy = [];
                                end
                                
                                norm_magn =  abs(SignalWithR2);
                                phase = angle(SignalWithR2);
                                norm_phase = phase - phase(1) - (phase(2)-phase(1))*(time - time(1))/(time(2) - time(1));
          
                                if dict_options.coordinate.classic_polar
                                    tempSignalVector.classic_polar((rot-1)*length_one_orientation + 1 : rot*length_one_orientation) = [theta_noisy norm_magn norm_phase];
                                end
                                
                                if dict_options.coordinate.classic_cartesian
                                    norm_complex_classic = norm_magn.*exp(1i*norm_phase);
                                    real_complex_classic = real(norm_complex_classic);
                                    imag_complex_classic = imag(norm_complex_classic);
                                    tempSignalVector.classic_cartesian((rot-1)*length_one_orientation + 1 : rot*length_one_orientation) = [theta_noisy real_complex_classic imag_complex_classic];
                                end
                                
                                if (dict_options.coordinate.polyfit_polar || dict_options.coordinate.polyfit_cartesian || dict_options.coordinate.polyfit_cartesian_demean)
                                    poly_coeff = polyfit(time, phase, 1);
                                    norm_phase_polyfit = phase - (time*poly_coeff(1) + poly_coeff(2));
   
                                    if dict_options.coordinate.polyfit_polar 
                                         tempSignalVector.polyfit_polar((rot-1)*length_one_orientation + 1 : rot*length_one_orientation) = [theta_noisy  norm_magn norm_phase_polyfit];
                                    end
                                    if dict_options.coordinate.polyfit_cartesian
                                        
                                        norm_complex_polyfit= norm_magn.*exp(1i*norm_phase_polyfit);
                                        real_complex_polyfit = real(norm_complex_polyfit);
                                        imag_complex_polyfit = imag(norm_complex_polyfit);
                                        
                                        tempSignalVector.polyfit_cartesian((rot-1)*length_one_orientation + 1 : rot*length_one_orientation) = [theta_noisy  real_complex_polyfit imag_complex_polyfit];
                                    end
                                    
                                    if dict_options.coordinate.polyfit_cartesian_demean
                                        norm_phase_polyfit_demean = norm_phase_polyfit - norm_phase_polyfit(1);
 
                                        norm_complex_polyfit_demean= norm_magn.*exp(1i*norm_phase_polyfit_demean);
                                        real_complex_polyfit_demean = real(norm_complex_polyfit_demean);
                                        imag_complex_polyfit_demean = imag(norm_complex_polyfit_demean);
                                        
                                        tempSignalVector.polyfit_cartesian_demean((rot-1)*length_one_orientation + 1 : rot*length_one_orientation) = [theta_noisy  real_complex_polyfit_demean imag_complex_polyfit_demean];
                                    end
                                end
                                
                            end

                            main_direction = SignalComponent(l,m,n,o,9).current_dir;

                            for kCoordinate = 1:lCoordinate
                                coordinate = coordinate_list{kCoordinate};
                                SignalValuesMultiCoordinates.(coordinate)(:,l,m,n,o,p,q,r) = tempSignalVector.(coordinate);
                            end
                            
                            gRatioDicoTemp(l,m,n,o,p,q,r) = Signal0.gRatio;
                            xiDicoTemp(l,m,n,o,p,q,r) = xi;
                            xaDicoTemp(l,m,n,o,p,q,r) = xa;
                            T2myelDicoTemp(l,m,n,o,p,q,r) = T2myel;
                            T2outDicoTemp(l,m,n,o,p,q,r) = T2out;
                            weightDicoTemp(l,m,n,o,p,q,r) = weight;
                            directionsDicoTemp(:,l,m,n,o,p,q,r) = main_direction;
                            thetaDicoTemp(l,m,n,o,p,q,r) = acos(abs(main_direction(3)));
                            
                        end
                    end
                end
            end
            toc()
        end
    end
end
toc()

gRatioValues = squeeze(gRatioDicoTemp);
xiValues = squeeze(xiDicoTemp);
xaValues = squeeze(xaDicoTemp);
T2myelValues = squeeze(T2myelDicoTemp);
T2outValues = squeeze(T2outDicoTemp);
weightValues = squeeze(weightDicoTemp);
directionsValues = squeeze(directionsDicoTemp);
thetaValues = squeeze(thetaDicoTemp);

gRatio = gRatioRange/100;

infoDico = 'In order, gRatio, xi, dir, T2myel, T2out, weight';
infoSignal = ['concatenation of ' num2str(nb_orientations) ' rotations each composed of the theta angle (B0 angle) 12 TE real and imag part of the signal with a polyfit normalization phase'];

if noise == 0.005
    prefix_name = 'SignalWithNoise05';
else
    prefix_name = ['SignalWithNoise' num2str(100*noise)];
end

if dict_options.include_theta
    suffix_theta = 'with_theta';
else
    suffix_theta = 'without_theta';
end

for kCoordinate = 1:lCoordinate
    coordinate = coordinate_list{kCoordinate};

    base_name = [prefix_name '_FVF' num2str(FVF) '_replic' num2str(replic) '_' num2str(nb_orientations) ...
                 'orientations_' length(time)  'TE_Porcine1_fix_xa_' coordinate '_'  suffix_theta];
    signal_name = [base_name '.h5py'];

    SignalValues = single(squeeze(SignalValuesMultiCoordinates.(coordinate)));
    
    save([output_folder '/' signal_name], 'SignalValues', ...
     'gRatioRange', 'xiRange', 'xaRange', 'T2myelRange','T2outRange','weightRange', ...
    'gRatioValues', 'xiValues', 'xaValues', 'T2myelValues', 'T2outValues', 'directionsValues', 'thetaValues', ...
    'weightValues', 'time', 'infoDico', 'infoSignal', 'sphere_rotations', '-v7.3')
end












