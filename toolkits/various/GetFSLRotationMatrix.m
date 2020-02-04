%% rot = GetFSLRotationMatrix(fn,path)
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 July 2018
% Date last modified:
%
%
function rot_mat = GetFSLRotationMatrix(fn,path)

if nargin==2
    fn = [path filesep fn];
end

cmd = ['avscale ' fn];

[~,cmdout] = system(cmd);

ind_begin   = strfind(cmdout,':');
ind_begin   = ind_begin(1)+2;
ind_end     = strfind(cmdout,'Scale');
ind_end     = ind_end(1) -1;

aff_mat = str2num(cmdout(ind_begin:ind_end));

rot_mat = aff_mat(1:3,1:3);

end