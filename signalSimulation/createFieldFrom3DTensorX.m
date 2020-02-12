function TotalField = BFieldFromTensorX3d(XTensor3d,B0,gamma,current_dir)
% keyboard;
dims = size(XTensor3d);

[kx,ky,kz] = ndgrid(-dims(1)/2:dims(1)/2-1, -dims(2)/2:dims(2)/2-1, -dims(3)/2:dims(3)/2-1);

kx = (kx / max(abs(kx(:)))) / dims(1);
ky = (ky / max(abs(ky(:)))) / dims(2);
kz = (kz / max(abs(kz(:)))) / dims(3);

kx = single(fftshift(kx));
ky = single(fftshift(ky));
kz = single(fftshift(kz));

k2 = kx.^2 + ky.^2 + kz.^2;
H_Vec = current_dir;
kH_over_k2 = (H_Vec(1) * kx + H_Vec(2) * ky + H_Vec(3) * kz) ./ (eps + k2);
clear k2

Fx(:,:,:,1) = fftn(XTensor3d(:,:,:,1));
Res =       ((H_Vec(1)^2)/3 - H_Vec(1)*kx .* kH_over_k2) .* Fx(:,:,:,1) ;
clear Fx
Fx(:,:,:,1) = fftn(XTensor3d(:,:,:,2));
Res =       Res+(2*(H_Vec(1)*H_Vec(2))/3 - (H_Vec(1)*ky + H_Vec(2)*kx) .* kH_over_k2) .* Fx(:,:,:,1) ;
clear Fx
Fx(:,:,:,1) = fftn(XTensor3d(:,:,:,3));
Res =       Res+(2*(H_Vec(1)*H_Vec(3))/3 - (H_Vec(1)*kz + H_Vec(3)*kx) .* kH_over_k2) .* Fx(:,:,:,1) ;
clear Fx kx
Fx(:,:,:,1) = fftn(XTensor3d(:,:,:,4));
Res =       Res+((H_Vec(2)^2)/3 - H_Vec(2)*ky .* kH_over_k2) .* Fx(:,:,:,1) ;
clear Fx
Fx(:,:,:,1) = fftn(XTensor3d(:,:,:,5));
Res =       Res+(2*(H_Vec(2)*H_Vec(3))/3 - (H_Vec(2)*kz + H_Vec(3)*ky) .* kH_over_k2) .* Fx(:,:,:,1) ;
clear Fx ky
Fx(:,:,:,1) = fftn(XTensor3d(:,:,:,6));
Res =       Res+((H_Vec(3)^2)/3 - H_Vec(3)*kz .* kH_over_k2) .* Fx(:,:,:,1);
clear Fx XTensor3d kz k2 kH_over_k2 


TotalField = ifftn(Res);
clear Res
TotalField = TotalField*gamma*B0;

end

%% OLD

%Reshape and fft XTensor
% for j = 1:6
%     Fx(:,:,:,j) = fftn(XTensor3d(:,:,:,j));
% end
% Fx(:,:,:,1) = fftn(XTensor3d(:,:,:,1));
% Res =       ((H_Vec(1)^2)/3 - H_Vec(1)*kx .* kH_over_k2) .* Fx(:,:,:,1) ;
% clear Fx
% Fx(:,:,:,2) = fftn(XTensor3d(:,:,:,2));
% Res =       Res+(2*(H_Vec(1)*H_Vec(2))/3 - (H_Vec(1)*ky + H_Vec(2)*kx) .* kH_over_k2) .* Fx(:,:,:,2) ;
% clear Fx
% Fx(:,:,:,3) = fftn(XTensor3d(:,:,:,3));
% Res =       Res+(2*(H_Vec(1)*H_Vec(3))/3 - (H_Vec(1)*kz + H_Vec(3)*kx) .* kH_over_k2) .* Fx(:,:,:,3) ;
% clear Fx
% Fx(:,:,:,4) = fftn(XTensor3d(:,:,:,4));
% Res =       Res+((H_Vec(2)^2)/3 - H_Vec(2)*ky .* kH_over_k2) .* Fx(:,:,:,4) ;
% clear Fx
% Fx(:,:,:,5) = fftn(XTensor3d(:,:,:,5));
% Res =       Res+(2*(H_Vec(2)*H_Vec(3))/3 - (H_Vec(2)*kz + H_Vec(3)*ky) .* kH_over_k2) .* Fx(:,:,:,5) ;
% clear Fx
% Fx(:,:,:,6) = fftn(XTensor3d(:,:,:,6));
% Res =       Res+((H_Vec(3)^2)/3 - H_Vec(3)*kz .* kH_over_k2) .* Fx(:,:,:,6);
% clear Fx XTensor3d


        
% Res =       ((H_Vec(1)^2)/3 - H_Vec(1)*kx .* kH_over_k2) .* Fx(:,:,:,1) + ...                         %   Fx11 
%             (2*(H_Vec(1)*H_Vec(2))/3 - (H_Vec(1)*ky + H_Vec(2)*kx) .* kH_over_k2) .* Fx(:,:,:,2) + ...    %   Fx12
%             (2*(H_Vec(1)*H_Vec(3))/3 - (H_Vec(1)*kz + H_Vec(3)*kx) .* kH_over_k2) .* Fx(:,:,:,3) + ...    %   Fx13
%             ((H_Vec(2)^2)/3 - H_Vec(2)*ky .* kH_over_k2) .* Fx(:,:,:,4) + ...                                    %   Fx22
%             (2*(H_Vec(2)*H_Vec(3))/3 - (H_Vec(2)*kz + H_Vec(3)*ky) .* kH_over_k2) .* Fx(:,:,:,5) + ...    %   Fx23
%             ((H_Vec(3)^2)/3 - H_Vec(3)*kz .* kH_over_k2) .* Fx(:,:,:,6);    
%         
