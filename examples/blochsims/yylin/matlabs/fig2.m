% ====== data for turbulent spin motion
format compact
% ==== const.
i=sqrt(-1);
gammah=2.67515255e8;
vacuum=12.566370614e-7;
hbar=1.05457266e-34;
avogd=6.0221367e23;
boltz=1.380658e-23;
epsx=1e-3;
% ===== physical para.
tr=10e-3;
conc=110*0.95;
% =================
larmor=600e6; %Hz
temp=300;
M0=gammah*hbar*tanh(hbar*pi*larmor/boltz/temp)*conc*avogd*1e3/2
td=1/(gammah*vacuum*M0) %sec
% ===== pulse
theta=90/180*pi; %pulse angle, rad
phi=90/180*pi; %pulse phase, rad
strg=5; %G/cm
durg=1e-3;  %sec
% ===== sample
nz=200;
nhelix=2;
ncypt=nz/nhelix; % # of pts in 1 helix cycle
ph=reshape((2*pi/100*[0:99].')*ones(1,nhelix),nz,1);
gammagt=durg*strg*gammah*1e-4;
lz=2*pi*nhelix/gammagt;
z=([0.5:-1/(nz-1):-0.5].')*lz;
% ====== equ. mag. 
M=zeros(nz,3);
size(M)
M(:,3)=ones(nz,1);
% ===== delta RF         
if theta ~= 0
   sthe=sin(theta);
   cthe=cos(theta);
   sphi=sin(phi);
   cphi=cos(phi);
   trans=zeros(3,3);
   trans(1,1)=cphi^2+sphi^2*cthe;
   trans(2,1)=sphi*cphi*(1-cthe);
   trans(3,1)=sphi*sthe;
   trans(1,2)=trans(2,1);
   trans(2,2)=sphi^2+cphi^2*cthe;
   trans(3,2)=-cphi*sthe;
   trans(1,3)=-trans(3,1);
   trans(2,3)=-trans(3,2);
   trans(3,3)=cthe;
   M=M*trans;
end
% ===== delta gradient
phstart=0.0d0/180*pi; %assume pulse phase 90 deg
sph=sin(ph+phstart);
cph=cos(ph+phstart);
m1=cph.*M(:,1)+sph.*M(:,2);
m2=-sph.*M(:,1)+cph.*M(:,2);
M(:,1)=m1;
M(:,2)=m2;
% ===== residual mag.
M(:,1)=M(:,1)+epsx;
% ===== free evolution
t=[0:0.01:1];
prec=1e-5;
[fidt,fid,Y]=fig2gen(t,M,gammah,td,tr,z,prec);
Y=reshape(Y,length(t),nz,3);

