function flirtAB2C(A,B,C,extension,varargin);
% function flirtAB2C(A,B,C,extension);
% A is input image where the registration will be done=A;
% B is the image to who the transformation is appliedOtherImage=B;
% C is the reference image reference=C;
% =varargin{1}; -cost normmi -dof 6 

input=A;
OtherImage=B;
reference=C;

if nargin==5
    opts=varargin{1};
else
    opts=[];
end
setenv('PATH', [getenv('PATH') ':/usr/local/fsl/bin']);

unix(['flirt ',opts,' -in ',input,' -ref ',reference,' -out ',input,extension,' -omat ',extension,'.mat']) 

if ~iscell(OtherImage)
    unix(['flirt -in ',OtherImage,' -ref ',reference,' -applyxfm -init ',extension,'.mat',' -out ',OtherImage,extension,' ']) 
else
    for k = 1:length(OtherImage)
        currentOtherImage = OtherImage{k};
        unix(['flirt -in ', currentOtherImage,' -ref ',reference,' -applyxfm -init ',extension,'.mat', ' -out ',currentOtherImage, extension,' ']) 
    end
end

% unix(['rm ',extension,'.mat'])

