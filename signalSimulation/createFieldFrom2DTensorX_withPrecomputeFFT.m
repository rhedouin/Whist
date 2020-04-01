function TotalField = createFieldFromTensorX2D_withPrecomputeFFT(kx,ky,k2,Fx,B0,H_Vec)
% From ApplySTI in COSMOS_STI_Toolbox
%%%% Short version
gamma = 42.6;

kH_over_k2 = (H_Vec(1) * kx + H_Vec(2) * ky) ./ (eps + k2);

Res =       ((H_Vec(1)^2)/3 - H_Vec(1)*kx .* kH_over_k2) .* Fx(:,:,1) + ...                         %   Fx11
    (2*(H_Vec(1)*H_Vec(2))/3 - (H_Vec(1)*ky + H_Vec(2)*kx) .* kH_over_k2) .* Fx(:,:,2) + ...    %   Fx12
    ((H_Vec(2)^2)/3 - H_Vec(2)*ky .* kH_over_k2) .* Fx(:,:,4) + ...                                 %   Fx22
    (H_Vec(3)^2)/3.* Fx(:,:,6);

TotalField = ifftn(Res);
TotalField = TotalField*gamma*B0;

end