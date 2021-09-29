% clear
close all
cd /project/3015069.04/temp/Jobs

req_mem   = 64e9;
req_etime = 20000;

fa_list = [10, 20, 30, 40, 50, 70];
for fa = fa_list
%     job = qsubfeval(@transformAndConcatenateBS3Data, fa, 'memreq',  req_mem,  'timreq',  req_etime);
    transformAndConcatenateBS3Data(fa);

end
