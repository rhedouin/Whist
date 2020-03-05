function out = computeAndSaveSignalFromModels(FVF_round , suffix, dico_inputfolder, dico_outputfolder, dict_params)

mkdir(dico_outputfolder);

gamma  = 42.6;

FVF_inputfolder = [dico_inputfolder 'FVF' num2str(FVF_round) '_N400' suffix '/'];
FVF_outputfolder = [dico_outputfolder 'FVF' num2str(FVF_round) '_N400' suffix '/'];
mkdir(FVF_outputfolder)

xiMyelinRange = dict_params.myelin.xiRange;
lXiMyelin = length(xiMyelinRange);

xaMyelinRange = dict_params.myelin.xaRange;
lXaMyelin = length(xaMyelinRange);

gRatioRange = dict_params.gRatioRange;
lGRatio = length(gRatioRange);

fiber_directions = dict_params.fiber_directions;
nb_fiber_directions = size(fiber_directions, 1);

TE = dict_params.TE;

if dict_params.rotations == 0
    sphere_rotations = [1 0 0; 0 1 0; 0 0 1];
    lRot = 1;
    display('no rotation')
else
    sphere_rotations = dict_params.sphere_rotations;
    lRot = size(sphere_rotations, 3);
    display('rotation')
end

if dict_params.dispersion == 0
    dispersion = 0;
    lDisp = 1;
    display('no dispersion')
else
    if ~isfield(dict_params, 'nb_orientation_for_dispersion')
        dict_params.nb_orientation_for_dispersion = 10;
    end
    
    dispersion_list = dict_params.dispersion_list;
    kappa_list = dict_params.kappa_list;
    lDisp = length(dispersion_list);
    display('dispersion')
end

display(['lGRatio : ' num2str(lGRatio)])
display(['lXiMyelin : ' num2str(lXiMyelin)])
display(['lXaMyelin : ' num2str(lXaMyelin)])
display(['nb_fiber_directions : ' num2str(nb_fiber_directions)])
display(['lRot : ' num2str(lRot)])
display(['lDisp : ' num2str(lDisp)])

gRatioValues = zeros(lGRatio, lXiMyelin, lXaMyelin, nb_fiber_directions, lRot, lDisp);
xiMyelinValues = zeros(lGRatio, lXiMyelin, lXaMyelin, nb_fiber_directions, lRot, lDisp);
xaMyelinValues = zeros(lGRatio, lXiMyelin, lXaMyelin, nb_fiber_directions, lRot, lDisp);

dispersionValues = zeros(lGRatio, lXiMyelin, lXaMyelin, nb_fiber_directions, lRot, lDisp);
directionValues = zeros(3, lGRatio, lXiMyelin, lXaMyelin, nb_fiber_directions, lRot, lDisp);

model_parameters.TE = dict_params.TE;
model_parameters.dispersion = dict_params.dispersion;
         
for k = 1: lGRatio
    k
    clear Model ZoomedModel axon_collection  
    gRatio = gRatioRange(k);

    axonFileName = [FVF_inputfolder 'FVF' num2str(FVF_round)  '_gRatio' num2str(gRatio) '_N400' suffix '.mat']
    load(axonFileName);
    
    model_parameters.FVF = FVF;
    model_parameters.g_ratio = g_ratio;
    model_parameters.dims = dims;
    model_parameters.mask = mask;
    
    for l = 1:lXiMyelin
        xiMyelin = xiMyelinRange(l);
        for m = 1:lXaMyelin
            xaMyelin = xaMyelinRange(m);
            
            ['gRatio ' num2str(k) ', xiMyelin : ' num2str(l) ', xaMyelin : ' num2str(m) ]
            

            model_parameters.myelin.xi = xiMyelin;
            model_parameters.myelin.xa = xaMyelin;

            [tensor_X, model]  = create2DTensorXFromAxonList(axon_collection,model_parameters);

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
                model_parameters.current_dir = main_dir;

                if and(dict_params.rotations == 0, dict_params.dispersion == 0)

                    field_complex  = createFieldFrom2DTensorX_withPrecomputeFFT(kx, ky, k2, Fx, dict_params.B0, main_dir);
                    field = real(field_complex);

                    signal_components(k, l, m, n, 1, 1) = reconstructSignalComponents(field, model, model_parameters);
                    directionValues(:, k, l, m, n, 1, 1) = main_dir;
                    
                elseif and(dict_params.rotations == 0, dict_params.dispersion == 1)                   

                    for o = 1:length(dispersion_list)

                        dispersion = dispersion_list(o);
                        kappa = kappa_list(o);
                        
                        if dispersion == 0
                            field_complex  = createFieldFrom2DTensorX_withPrecomputeFFT(kx, ky, k2, Fx, dict_params.B0, main_dir);
                        else
                            [orientations, dispersion] = generateVMFsample_for_dispersion(main_dir, kappa, dispersion, dict_params.nb_orientation_for_dispersion);
                            
                            for p = 1:length(orientations)
                                field_complex(:,:,p)  = createFieldFromTensorX2D_withPrecomputeFFT(kx, ky, k2, Fx, dict_params.B0, orientations(p,:));
                            end
                        end
                        
                        field = real(field_complex);          
                        signal_components(k, l, m, n, 1, o) = reconstructSignalComponents(field, model, model_parameters);

                        clear field field_complex
                        directionValues(:, k, l, m, n, 1, o) = main_dir;

                    end
 
                elseif and(dict_params.rotations == 1, dict_params.dispersion == 0)

                    for o = 1:lRot
                                                
                        current_dir =  sphere_rotations(:,:,o) * main_dir;
                        model_parameters.current_dir = current_dir;

                        field_complex  = createFieldFrom2DTensorX_withPrecomputeFFT(kx, ky, k2, Fx, dict_params.B0, current_dir);
                        field = real(field_complex);
                        
                        signal_components(k, l, m, n, o, 1) = reconstructSignalComponents(field, model, model_parameters);
                        directionValues(:, k, l, m, n, o, 1) = current_dir;
                    end
                    
                elseif and(dict_params.rotations == 1, dict_params.dispersion == 1)

                    for o = 1:lRot                        
                        current_dir =  sphere_rotations(:,:,o) * main_dir;
                        model_parameters.current_dir = current_dir;
                   
                        for p = 1:length(dispersion_list)
                           ['nb_dir ' num2str(n) ', nb_rot ' num2str(o) ', nb_disp ' num2str(p)]

                            dispersion = dispersion_list(p);
                            kappa = kappa_list(p);
                            
                            if dispersion == 0
                                field_complex  = createFieldFrom2DTensorX_withPrecomputeFFT(kx, ky, k2, Fx, dict_params.B0, main_dir);
                            else
                                [orientations, dispersion] = generateVMFsample_for_dispersion(main_dir, kappa, dispersion, dict_params.nb_orientation_for_dispersion);
                                
                                for q = 1:length(orientations)
                                    field_complex(:,:,q)  = createFieldFrom2DTensorX_withPrecomputeFFT(kx, ky, k2, Fx, dict_params.B0, orientations(p,:));
                                end
                            end
                            
                            field = real(field_complex);
                            clear field_complex

                            signal_components(k, l, m, n, o, p) = reconstructSignalComponents(field, model, model_parameters);
                            directionValues(:, k, l, m, n, o, p) = current_dir;
                            
                        end
                    end
                end
            end
            
            gRatioValues(k,l,m,n,:,:) = gRatio;
            xiMyelinValues(k,l,m,n,:,:) = xiMyelin;
            xaMyelinValues(k,l,m,n,:,:) = xaMyelin;     
            dispersionValues(k,l,m,n,:,:) = dispersion;

        end
    end
end

thetaValues = acos(abs(directionValues(3, :, :, :, :, :, :)));
info = 'In order, gRatio, xi, xa, fiber directions, rotation, dispersion';

save([FVF_outputfolder 'Signal_FVF' num2str(FVF_round) suffix '.mat'], 'signal_components', 'gRatioRange','gRatioValues', 'xiMyelinRange',  'xiMyelinValues', 'xaMyelinRange',  'xaMyelinValues', 'directionValues', 'thetaValues', 'dispersionValues', 'TE', 'info', 'sphere_rotations', 'fiber_directions', 'dims', 'dict_params');

out = 0
end

