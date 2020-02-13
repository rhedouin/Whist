function TotalField = createFieldFrom2DTensorX(tensor_X, model_parameters)

gamma = 42.6;

dims = size(tensor_X);
[kx,ky] = ndgrid(-dims(1)/2:dims(1)/2-1, -dims(2)/2:dims(2)/2-1);

kx = (kx / max(abs(kx(:)))) / dims(1);
ky = (ky / max(abs(ky(:)))) / dims(2);

kx = fftshift(kx);
ky = fftshift(ky);

k2 = kx.^2 + ky.^2;

%Reshape and fft XTensor
for j = 1:6
    Fx(:,:,j) = fftn(tensor_X(:,:,j));
end

H_Vec = model_parameters.field_direction;

kH_over_k2 = (H_Vec(1) * kx + H_Vec(2) * ky) ./ (eps + k2);

Res =       ((H_Vec(1)^2)/3 - H_Vec(1)*kx .* kH_over_k2) .* Fx(:,:,1) + ...                         %   Fx11 
    (2*(H_Vec(1)*H_Vec(2))/3 - (H_Vec(1)*ky + H_Vec(2)*kx) .* kH_over_k2) .* Fx(:,:,2) + ...    %   Fx12
    (2*(H_Vec(1)*H_Vec(3))/3 - H_Vec(3)*kx .* kH_over_k2) .* Fx(:,:,3) + ...    %   Fx13
    ((H_Vec(2)^2)/3 - H_Vec(2)*ky .* kH_over_k2) .* Fx(:,:,4) + ...                                    %   Fx22
    (2*(H_Vec(2)*H_Vec(3))/3 - H_Vec(3)*ky.* kH_over_k2) .* Fx(:,:,5) + ...    %   Fx23
    (H_Vec(3)^2)/3.* Fx(:,:,6);        %   Fx33

TotalField = ifftn(Res);
TotalField = TotalField*gamma*model_parameters.B0;
end