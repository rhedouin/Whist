function [ratio, weight_out, weight_myel] = computeRelativeWeightFromT1Effect(T1_intra_axonal, T1myel, rho, TR, flip_angle)

E1_intra_axonal = exp(-TR/T1_intra_axonal);
E1myel = exp(-TR/T1myel);

weight_out = sind(flip_angle) * (1 - E1out) / (1 - cosd(flip_angle) * E1out);
weight_out = sind(flip_angle) * (1 - E1out) / (1 - cosd(flip_angle) * E1out);

weight_myel = sind(flip_angle) * (1 - E1myel) / (1 - cosd(flip_angle) * E1myel);

ratio =  (rho*weight_myel) / weight_out;
end