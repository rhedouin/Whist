function out = createDictionaryPartSignalMultiModalities(signal_path, output_folder, experience_name, T2MyelinRange, T2IntraExtraAxonalRange, weightRange, nb_TE, noise, FVF, nb_orientations, replic, dict_options)
load(signal_path); 
TE = TE(1:nb_TE); 

lgRatio = length(gRatioRange); 
lXiMyelin = length(xiMyelinRange); 
lXaMyelin = length(xaMyelinRange); 
lDir = size(fiber_directions,1);
lT2Myelin = length(T2MyelinRange);
lT2IntraExtraAxonal = length(T2IntraExtraAxonalRange);
lWeight = length(weightRange); 

if dict_options.include_theta
   length_one_orientation = 2*nb_TE + 1;
else
   length_one_orientation = 2*nb_TE;
end
lVec = nb_orientations * length_one_orientation;

coordinate_total_list = fieldnames(dict_options.coordinate);
it = 0;

for kCoordinate = 1:length(coordinate_total_list)
    if dict_options.coordinate.(coordinate_total_list{kCoordinate}) == 1
        it = it + 1;
        coordinate_list{it} = coordinate_total_list{kCoordinate};
        SignalValuesMultiCoordinates.(coordinate_list{it}) = zeros(lVec, lgRatio, lXiMyelin, lXaMyelin, lDir, lT2Myelin, lT2IntraExtraAxonal, lWeight, 'single');
    end
end
lCoordinate = length(coordinate_list);

gRatioValues = zeros(lgRatio, lXiMyelin, lXaMyelin, lDir, lT2Myelin, lT2IntraExtraAxonal, lWeight, 'single');
xiMyelinValues = zeros(lgRatio, lXiMyelin, lXaMyelin, lDir, lT2Myelin, lT2IntraExtraAxonal, lWeight, 'single');
xaMyelinValues = zeros(lgRatio, lXiMyelin, lXaMyelin, lDir, lT2Myelin, lT2IntraExtraAxonal, lWeight, 'single');
T2MyelinValues = zeros(lgRatio, lXiMyelin, lXaMyelin, lDir, lT2Myelin, lT2IntraExtraAxonal, lWeight, 'single');
T2IntraExtraAxonalValues = zeros(lgRatio, lXiMyelin, lXaMyelin, lDir, lT2Myelin, lT2IntraExtraAxonal, lWeight, 'single');
weightValues = zeros(lgRatio, lXiMyelin, lXaMyelin, lDir, lT2Myelin, lT2IntraExtraAxonal, lWeight, 'single');
directionsValues = zeros(3,lgRatio, lXiMyelin, lXaMyelin, lDir, lT2Myelin, lT2IntraExtraAxonal, lWeight, 'single');
thetaValues = zeros(lgRatio, lXiMyelin, lXaMyelin, lDir, lT2Myelin, lT2IntraExtraAxonal, lWeight, 'single');

tic()
for l = 1:lgRatio
    gRatio = gRatioRange(l)
    for m = 1:lXiMyelin
        xiMyelin = xiMyelinRange(m)
        for n = 1:lXaMyelin
            xaMyelin = xaMyelinRange(n)
            tic()
            for o = 1:lDir
                for p = 1:lT2Myelin
                    T2Myelin = T2MyelinRange(p);
                    for q = 1:lT2IntraExtraAxonal
                        T2IntraExtraAxonal = T2IntraExtraAxonalRange(q);
                        for r = 1:lWeight
                            weight = weightRange(r);
                            
                            for kCoordinate = 1:lCoordinate
                                coordinate = coordinate_list{kCoordinate};
                                tempSignalVector.(coordinate) = zeros(lVec, 1);
                            end
                                                        
                            for rot = 1:nb_orientations
                                Signal0 = signal_components(l,m,n,o,rot);

                                SignalWithR2 = exp(-TE/T2IntraExtraAxonal).*Signal0.intra_axonal.signal(1:nb_TE) + weight*exp(-TE/T2Myelin).*Signal0.myelin.signal(1:nb_TE) + ...
                                               exp(-TE/T2IntraExtraAxonal).*Signal0.extra_axonal.signal(1:nb_TE);
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
                                norm_phase = phase - phase(1) - (phase(2)-phase(1))*(TE - TE(1))/(TE(2) - TE(1));
          
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
                                    poly_coeff = polyfit(TE, phase, 1);
                                    norm_phase_polyfit = phase - (TE*poly_coeff(1) + poly_coeff(2));
   
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

                            main_direction = signal_components(l,m,n,o,6).current_dir;

                            for kCoordinate = 1:lCoordinate
                                coordinate = coordinate_list{kCoordinate};
                                SignalValuesMultiCoordinates.(coordinate)(:,l,m,n,o,p,q,r) = tempSignalVector.(coordinate);
                            end
                            
                            gRatioValues(l,m,n,o,p,q,r) = Signal0.g_ratio;
                            xiMyelinValues(l,m,n,o,p,q,r) = xiMyelin;
                            xaMyelinValues(l,m,n,o,p,q,r) = xaMyelin;
                            T2MyelinValues(l,m,n,o,p,q,r) = T2Myelin;
                            T2IntraExtraAxonalValues(l,m,n,o,p,q,r) = T2IntraExtraAxonal;
                            weightValues(l,m,n,o,p,q,r) = weight;
                            directionsValues(:,l,m,n,o,p,q,r) = main_direction;
                            thetaValues(l,m,n,o,p,q,r) = acos(abs(main_direction(3)));
                            
                        end
                    end
                end
            end
            toc()
        end
    end
end
toc()

gRatio = gRatioRange/100;

infoDico = 'In order, gRatio, xiMyelin, xaMyelin, dir, T2Myelin, T2IntraExtraAxonal, weight';
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
                 'orientations_' num2str(nb_TE)  'TE_' experience_name '_fix_xa_' coordinate '_'  suffix_theta];
    signal_name = [base_name '.h5py'];

    SignalValues = single(SignalValuesMultiCoordinates.(coordinate));
    
    save([output_folder '/' signal_name], 'SignalValues', ...
     'gRatioRange', 'xiMyelinRange', 'xaMyelinRange', 'T2MyelinRange','T2IntraExtraAxonalRange','weightRange', ...
    'gRatioValues', 'xiMyelinValues', 'xaMyelinValues', 'T2MyelinValues', 'T2IntraExtraAxonalValues', 'directionsValues', 'thetaValues', ...
    'weightValues', 'TE', 'infoDico', 'infoSignal', 'sphere_rotations', '-v7.3')

end
out = 0;
end











