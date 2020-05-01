function out = createRotationMatrixWithArbitraryCenter(R, c)
c_neg = [1 0 0 -c(1); 0 1 0 -c(2); 0 0 1 -c(3); 0 0 0 1];
c_pos = [1 0 0 c(1); 0 1 0 c(2); 0 0 1 c(3); 0 0 0 1];

R0 = zeros(4);
R0(1:3, 1:3) = R;
R0(4, 4) = 1;

out = c_pos*R0*c_neg;

end