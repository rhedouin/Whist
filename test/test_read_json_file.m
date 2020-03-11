% test read json file

folder = '/project/3015069.04/deep_learning/multi_orientations/Porcine2/lowres/SignalWithNoise05_5rep_6orientations_18TE_Porcine2_fix_xa_large_FVF_20_directions_polyfit_cartesian_with_theta_with0reg/';
json_file = [folder 'total_history.json'];

fid = fopen(json_file);
raw = fread(fid,inf);
str = char(raw'); 
