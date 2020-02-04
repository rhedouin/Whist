function out = computeAndSaveSignalFromModels(FVF , suffix, dico_inputfolder, dico_outputfolder, options)

mkdir(dico_outputfolder);

B0 = options.B0;
gamma  = 42.6;

FVF_inputfolder = [dico_inputfolder 'FVF' num2str(FVF) '_N400' suffix '/'];
FVF_outputfolder = [dico_outputfolder 'FVF' num2str(FVF) '_N400' suffix '/'];
mkdir(FVF_outputfolder)

time = options.time;

xiRange = options.xiRange;
lXi = length(xiRange);

xaRange = options.xaRange;
lXa = length(xaRange);

gRatioRange = options.gRatioRange;
lGRatio = length(gRatioRange);

mask = zeros(1000, 1000);
mask(round(1000/3):round(2*1000/3), round(1000/3):round(2*1000/3)) = 1;

fiber_directions = options.fiber_directions;
nb_fiber_directions = size(fiber_directions, 1);

if options.rotations == 0
    sphere_rotations = [1 0 0; 0 1 0; 0 0 1];
    lRot = 1;
    display('no rotation')
else
    sphere_rotations = options.sphere_rotations;
    lRot = size(sphere_rotations, 3);
    display('rotation')
end

if options.dispersion == 0
    dispersion = 0;
    lDisp = 1;
    display('no dispersion')
else
    if isfield(options, 'nb_orientation_for_dispersion')
        nb_orientation_for_dispersion = options.nb_orientation_for_dispersion;
    else
        nb_orientation_for_dispersion = 10;
    end
    
    dispersion_list = options.dispersion_list;
    kappa_list = options.kappa_list;
    lDisp = length(dispersion_list);
    display('dispersion')
end

display(['lGRatio : ' num2str(lGRatio)])
display(['lXi : ' num2str(lXi)])
display(['lXa : ' num2str(lXa)])
display(['nb_fiber_directions : ' num2str(nb_fiber_directions)])
display(['lRot : ' num2str(lRot)])
display(['lDisp : ' num2str(lDisp)])

SignalComponent = repmat(struct('Axon', 1, 'Myelin', 1, 'Extra', 1, 'FVF', 1, 'gRatio', 1, ...
    'xi', 1, 'xa', 1, 'current_dir', 1, 'dispersion', 1, 'time', 1), ...
    [lGRatio, lXi, lXa, nb_fiber_directions, lRot, lDisp]);

gRatioValues = zeros(lGRatio, lXi, lXa, nb_fiber_directions, lRot, lDisp);
xiValues = zeros(lGRatio, lXi, lXa, nb_fiber_directions, lRot, lDisp);
xaValues = zeros(lGRatio, lXi, lXa, nb_fiber_directions, lRot, lDisp);
dispersionValues = zeros(lGRatio, lXi, lXa, nb_fiber_directions, lRot, lDisp);
directionValues = zeros(3, lGRatio, lXi, lXa, nb_fiber_directions, lRot, lDisp);

for k = 1: lGRatio
    k
    clear Model ZoomedModel axon_collection  
    gRatio = gRatioRange(k);

    axonFileName = [FVF_inputfolder 'AxonMap_FVF' num2str(FVF)  '_gRatio' num2str(gRatio) '_N400' suffix '.mat']
    load(axonFileName);

    side = dims(1);
    
    for l = 1:lXi
        xi = xiRange(l);
        for m = 1:lXa
            xa = xaRange(m);
            
            ['gRatio ' num2str(k) ', xi ' num2str(l)]

            [tensor_X, model]  = create2DTensorXFromAxonList(axon_collection,dims,xa,xi);

            Fx(:,:,1) = fftn(tensor_X(:,:,1));
            Fx(:,:,2) = fftn(tensor_X(:,:,2));
            Fx(:,:,4) = fftn(tensor_X(:,:,4));
            Fx(:,:,6) = fftn(tensor_X(:,:,6));

            [kx,ky] = ndgrid(-dims(1)/2:dims(1)/2-1, -dims(2)/2:dims(2)/2-1);

            kx = (kx / max(abs(kx(:)))) / dims(1);
            ky = (ky / max(abs(ky(:)))) / dims(2);

            kx = fftshift(kx);
            ky = fftshift(ky);

            k2 = kx.^2 + ky.^2;
            
            for n = 1:nb_fiber_directions
                ['nb_dir ' num2str(n)]

                main_dir =  fiber_directions(n,:)';
                if and(options.rotations == 0, options.dispersion == 0)

                    field_complex  = createFieldFrom2DTensorX_withPrecomputeFFT(kx, ky, k2, Fx, B0, gamma, main_dir);
                    field = real(field_complex);

                    SignalComponent(k, l, m, n, 1, 1) = reconstructSignalComponents(field, model, mask, time, current_FVF, current_g_ratio, xi, xa, main_dir, dispersion);
                    directionValues(:, k, l, m, n, 1, 1) = main_dir;
                    
                elseif and(options.rotations == 0, options.dispersion == 1)                   

                    for o = 1:length(dispersion_list)
                        dispersion = dispersion_list(o);
                        kappa = kappa_list(o);
                        
                        if dispersion == 0
                            field_complex  = createFieldFrom2DTensorX_withPrecomputeFFT(kx, ky, k2, Fx, B0, gamma, main_dir);
                        else
                            [orientations, dispersion] = generateVMFsample_for_dispersion(main_dir, kappa, dispersion, nb_orientation_for_dispersion);
                            
                            for p = 1:length(orientations)
                                field_complex(:,:,p)  = createFieldFromTensorX2D_withPrecomputeFFT(kx, ky, k2, Fx, B0, gamma, orientations(p,:));
                            end
                        end
                        
                        field = real(field_complex);          
                        SignalComponent(k, l, m, n, 1, o) = reconstructSignalComponents(field, model, mask, time, current_FVF, current_g_ratio, xi, xa, main_dir, dispersion);

                        clear field field_complex
                        directionValues(:, k, l, m, n, 1, o) = main_dir;

                    end
 
                elseif and(options.rotations == 1, options.dispersion == 0)

                    for o = 1:lRot
                                                
                        current_dir =  sphere_rotations(:,:,o) * main_dir;
                        
                        field_complex  = createFieldFrom2DTensorX_withPrecomputeFFT(kx, ky, k2, Fx, B0, gamma, current_dir);
                        field = real(field_complex);
                        
                        SignalComponent(k, l, m, n, o, 1) = reconstructSignalComponents(field, model, mask, time, current_FVF, current_g_ratio, xi, xa, current_dir, dispersion);
                        directionValues(:, k, l, m, n, o, 1) = current_dir;

                    end
                    
                elseif and(options.rotations == 1, options.dispersion == 1)

                    for o = 1:lRot                        
                        current_dir =  sphere_rotations(:,:,o) * main_dir;
                        
                        for p = 1:length(dispersion_list)
                           ['nb_dir ' num2str(n) ', nb_rot ' num2str(o) ', nb_disp ' num2str(p)]

                            dispersion = dispersion_list(p);
                            kappa = kappa_list(p);
                            
                            if dispersion == 0
                                field_complex  = createFieldFrom2DTensorX_withPrecomputeFFT(kx, ky, k2, Fx, B0, gamma, main_dir);
                            else
                                [orientations, dispersion] = generateVMFsample_for_dispersion(main_dir, kappa, dispersion, nb_orientation_for_dispersion);
                                
                                for q = 1:length(orientations)
                                    field_complex(:,:,q)  = createFieldFrom2DTensorX_withPrecomputeFFT(kx, ky, k2, Fx, B0, gamma, orientations(p,:));
                                end
                            end
                            
                            field = real(field_complex);
                            clear field_complex

                            SignalComponent(k, l, m, n, o, p) = reconstructSignalComponents(field, model, mask, time, current_FVF, current_g_ratio, xi, xa, current_dir, dispersion);
                            directionValues(:, k, l, m, n, o, p) = current_dir;
                            
                        end
                    end
                end
            end
            
            gRatioValues(k,l,m,n,:,:) = gRatio;
            xiValues(k,l,m,n,:,:) = xi;
            xaValues(k,l,m,n,:,:) = xa;     
            dispersionValues(k,l,m,n,:,:) = dispersion;

        end
    end
end

thetaValues = acos(abs(directionValues(3, :, :, :, :, :, :)));

info = 'In order, gRatio, xi, xa, fiber directions, rotation, dispersion';

save([FVF_outputfolder 'Signal_FVF' num2str(FVF) suffix '.mat'], 'SignalComponent', 'gRatioRange','gRatioValues', 'xiRange',  'xiValues', 'xaRange',  'xaValues', 'directionValues', 'thetaValues', 'dispersionValues', 'time', 'info', 'sphere_rotations', 'fiber_directions', 'dims', 'options');

out = 0
end

