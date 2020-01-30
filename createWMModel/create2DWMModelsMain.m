% clear
close all
cd /project/3015069.04/Jobs
basefolder = '/project/3015069.01/model/data/dicos/';

req_mem   = 6e9;
req_memoverhead   = 6e9;
req_etime = 5000;

it = 0;
expected_g_ratio = 0.5 : 0.05 : 0.85;
for k = 5 : 8
    suffix = ['_train' num2str(k)];
    for expected_FVF = 0.1 : 0.1 : 0.2
        it = it +1;

        if k < 5
            mode = 'remove'
        else
            mode = 'spread'
        end
        createOne2DWMModel(expected_FVF, expected_g_ratio, suffix, mode);

%         job{it} = qsubfeval(@createSave2DWMModels, expected_FVF, expected_g_ratio, suffix, mode, 'memreq',  req_mem,  'timreq',  req_etime, 'memoverhead',req_memoverhead);
    end
end

