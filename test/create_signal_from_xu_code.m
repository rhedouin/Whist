%% Field perturbation calculation using the Fourier method (Liu 2010) %% and signal calculation
clear;
clc;
% close all 
%% variables
%angle between B0 and long axis of thenested cylinder 
theta = pi/2;
sus_my=-10*10^-8; %magnitude ofisotropic and anisotropic susceptibility 
t2my=0.015; %t2 ofmyelin 
t2axea=0.050; %t2 ofintra/extra axonal 
density=0.7; %density 
gratio=0.65; %g-ratio
L=sqrt((100^2).*density./pi); %space is100 by100.

L=L/10;
L1=gratio.*L; 
gamma=42.576*10^6; %Hz/Tesla 
B0=3; %fieldstrength 
CONSTANT=B0*gamma*sus_my; 

%% meshgrid
sc=0.03; [x y]=meshgrid(-50:sc:50); 
px=0; 
py=0; 
pz=0; 
IX=x; 
IY=y; 
IZ=1; 
%% masks

R= sqrt((x-px).^2+(y-py).^2); 
R(R>L)=0; 
R(R<L1)=0;
R(R~=0)=sus_my; 
MASK1=logical(R);

%% angle phi 
phi=atan2( (y-py),(x-px) ); 
phi=mod(-phi,2*pi); %this works

a=123;
phi=MASK1.*phi;
%% ISO prep the tensor 
AA=1; BB=1; CC=1; %ISOTROPIC 

TH=sin(theta);
PA=((AA*(cos(phi).^2))+(BB*(sin(phi).^2)))*sus_my; 
PB=((-AA*cos(phi).*sin(phi)) + (BB*cos(phi).*sin(phi)) )*sus_my;

PA=PA.*MASK1; 
PB=PB.*MASK1;

%componentone 
partONE=fftshift(fftn(fftshift( PA ))).*(1/3);

%componenttwo algorithm 
FTB=fftshift(fftn(fftshift(PB.*TH))); 
FTA=fftshift(fftn(fftshift(PA.*TH))); 
cy=0; 
cx=0;
kx=(IY-cy).*FTB; 
ky=(IX-cx).*FTA; 
tot=kx+ky;

partTWO= ((TH.*(IX-cx))./((IX-cx).^2 + (IY-cy).^2)) .* tot;

%combineandsolve 
IN=partONE-partTWO;
ISO=ifftshift(ifftn(ifftshift(IN)))*B0*gamma;

%% ANISO prep the tensor 
AA=1; BB=-0.5; %ISOTROPIC 
TH=sin(theta);
PA=((AA*(cos(phi).^2))+(BB*(sin(phi).^2)))*sus_my; 
PB=((-AA*cos(phi).*sin(phi)) + (BB*cos(phi).*sin(phi)) )*sus_my;

PA=PA.*MASK1; 
PB=PB.*MASK1;

%component one
partONE=fftshift(fftn(fftshift( PA - ((cos(theta)^2).*(PA))))).*(1/3);

%componenttwo algorithm 
FTB=fftshift(fftn(fftshift(PB.*TH))); 
FTA=fftshift(fftn(fftshift(PA.*TH))); 
cy=0; 
cx=0;
kx=(IY-cy).*FTB; 
ky=(IX-cx).*FTA; 
tot=kx+ky;
partTWO= ((TH.*(IX-cx) + cos(theta))./((IX-cx).^2 + (IY-cy).^2 )) .* tot;

%combine and solve 
IN=partONE-partTWO;
ANISO=ifftshift(ifftn(ifftshift(IN)))*B0*gamma; 

%%
close all;
RT=real(ISO)+real(ANISO); 

figure; 
imagesc(real(ISO))
title('xu ISO')
colorbar
figure; 
imagesc(real(ANISO))
title('xu ANISO')
colorbar

figure;imagesc(squeeze(RT(:,:)));colorbar;title('total field (hz)');

%% masks
R= sqrt((x-px).^2+(y-py).^2); 
Rax=R; 
R2=R; 
R1=R;
R(R>L)=0; 
R(R<L1)=0; 
R(R~=0)=1;
MASK=logical(R); 
MASK=single(MASK); 
MASK(MASK==0)=nan; 
MASK=double(MASK);

R1(R1>L)=0;

M1=logical(R1); 
MASK1=MASK;
MASK1(isnan(MASK1))=0; 
MASK1=double(MASK1);

clear MASK 
R2(R2<L)=1; 
R2(R2~=1)=0;

Rax(Rax<L1)=1; 
Rax(Rax~=1)=0; 
%% EA mask

Rea=zeros(3334,3334); 
Rea(1500:1500+333,1500:1500+333)=1; 
TT=R+Rax; 
Rea=Rea-TT; 
clear TT MASK1 M1 R1

%% SIGNAL EVO 
tlimit=0.08; 
samples=100;
%% myelin signal

time=linspace(0.0001,tlimit,samples); 
lin_my=nonzeros(RT.*R);

allSIGmy=complex(zeros(1,samples)); 
for t=1:numel(time) 
    %tic
    tim=time(t);
    signal=exp((2*pi*tim*1i).*lin_my); 
    sumsig=sum(signal);
    CS=sumsig*exp(-tim/t2my)*sc*sc; 
    CS=CS*0.5; %proton density 
    allSIGmy(:,t)=CS;
end

%% intra-axonal signal

time=linspace(0.0001,tlimit,samples); 
lin_ax=nonzeros(RT.*Rax); 
allSIGax=complex(zeros(1,samples)); 

for t=1:numel(time)

    tim=time(t);
    signal=exp((2*pi*tim*1i).*lin_ax); 
    sumsig=sum(signal);
    CS=sumsig*exp(-tim/t2axea)*sc*sc; 
    
    allSIGax(:,t)=CS;
end

%% extra-axonal signal 
time=linspace(0.0001,tlimit,samples);

lin_ea=nonzeros(RT.*Rea); 
allSIGEA=complex(zeros(1,samples)); 
for t=1:numel(time)

    tim=time(t);
    signal=exp((2*pi*tim*1i).*lin_ea); 
    sumsig=sum(signal);

    CS=sumsig*exp(-tim/t2axea)*sc*sc; 
    
    allSIGEA(:,t)=CS;
end

tot=allSIGEA+allSIGax+allSIGmy; 
phase=angle(tot); 
mag=abs(tot); 
%%
figure;
% subplot(121);plot(time,mag/mag(1));title('magnitude evolution') 
subplot(121);plot(time,mag);title('magnitude evolution') 
subplot(122);plot(time,unwrap(phase));title('phase evolution')
return;






