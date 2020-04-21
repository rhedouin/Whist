function TotalField = createFieldFrom3DTensorX(XTensor3d, model_parameters)

gamma = 42.6;

dims = size(XTensor3d);

[kx,ky,kz] = ndgrid(-dims(1)/2:dims(1)/2-1, -dims(2)/2:dims(2)/2-1, -dims(3)/2:dims(3)/2-1);

kx = (kx / max(abs(kx(:)))) / dims(1);
ky = (ky / max(abs(ky(:)))) / dims(2);
kz = (kz / max(abs(kz(:)))) / dims(3);

kx = single(fftshift(kx));
ky = single(fftshift(ky));
kz = single(fftshift(kz));

k2 = kx.^2 + ky.^2 + kz.^2;
H_Vec = model_parameters.field_direction;
kH_over_k2 = (H_Vec(1) * kx + H_Vec(2) * ky + H_Vec(3) * kz) ./ (eps + k2);
clear k2

display('fft1')
Fx = fftn(XTensor3d(:,:,:,1));
Res =       ((H_Vec(1)^2)/3 - H_Vec(1)*kx .* kH_over_k2) .* Fx ;
clear Fx

display('fft2')
Fx = fftn(XTensor3d(:,:,:,2));
Res =       Res+(2*(H_Vec(1)*H_Vec(2))/3 - (H_Vec(1)*ky + H_Vec(2)*kx) .* kH_over_k2) .* Fx;
clear Fx

display('fft3')
Fx= fftn(XTensor3d(:,:,:,3));
Res =       Res+(2*(H_Vec(1)*H_Vec(3))/3 - (H_Vec(1)*kz + H_Vec(3)*kx) .* kH_over_k2) .* Fx;
clear Fx kx

display('fft4')
Fx = fftn(XTensor3d(:,:,:,4));
Res =       Res+((H_Vec(2)^2)/3 - H_Vec(2)*ky .* kH_over_k2) .* Fx;
clear Fx

display('fft5')
Fx = fftn(XTensor3d(:,:,:,5));
Res =       Res+(2*(H_Vec(2)*H_Vec(3))/3 - (H_Vec(2)*kz + H_Vec(3)*ky) .* kH_over_k2) .* Fx;
clear Fx ky

display('fft6')
Fx = fftn(XTensor3d(:,:,:,6));
Res =       Res+((H_Vec(3)^2)/3 - H_Vec(3)*kz .* kH_over_k2) .* Fx;
clear Fx XTensor3d kz k2 kH_over_k2 

display('ifft')
TotalField = ifftn(Res);
clear Res
TotalField = TotalField*gamma*model_parameters.B0;

end
