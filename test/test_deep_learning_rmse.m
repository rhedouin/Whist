% label:  FVF , rmse:  0.107789237594066
% label:  gRatio , rmse:  0.05782324344298481
% label:  xi , rmse:  0.025840281871222434
% label:  T2out , rmse:  0.007078574894693423
% label:  T2myel , rmse:  0.001075436861704524
% label:  weight , rmse:  0.4751495355988897
% 
% label: FVF, rmse: 0.11071181687113087
% label: gRatio, rmse: 0.0625974033966443
% label: xi, rmse: 0.02679734357945411
% label: xa, rmse: 0.0033305297101334347
% label: T2out, rmse: 0.007062677949621865
% label: T2myel, rmse: 0.0014170736455215813
% label: weight, rmse: 0.48697835805912254

rmse_without_xa = [0.107789237594066, 0.05782324344298481, 0.025840281871222434, 0.007078574894693423, 0.001075436861704524, 0.4751495355988897];
rmse_with_xa = [0.11071181687113087, 0.0625974033966443, 0.02679734357945411, 0.0033305297101334347, 0.007062677949621865, 0.0014170736455215813, 0.48697835805912254];

output_scale_without_xa = [0.85, 0.8511078953742981, 0.20000000298023224, 0.10000000149011612, 0.019999999552965164, 3.0];
output_scale_with_xa = [0.85, 0.8511078953742981, 0.20000000298023224, 0.10000000149011612, 0.10000000149011612, 0.019999999552965164, 3.0];

norm_rmse_without_xa = rmse_without_xa ./ output_scale_without_xa
norm_rmse_with_xa = rmse_with_xa ./ output_scale_with_xa

mean_norm_rmse_without_xa = mean(norm_rmse_without_xa)
mean_norm_rmse_with_xa = mean(norm_rmse_with_xa)

range = [0.45, 0.35, 0.4, 0.8, 0.16, 2.5]

norm_rmse_without_xa_range = norm_rmse_without_xa ./ range
mean_norm_rmse_without_xa_range = mean(norm_rmse_without_xa_range)
