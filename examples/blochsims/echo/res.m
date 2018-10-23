

w=25*pi*2;
tau=0.082;
Mo=1e-6;
phi=0*pi/180;
theta=180*pi/180;
D=1e8;

npts=256;
t2=(0:npts)/npts;
npts2d=256;
t1=tau:tau:(npts2d*tau);
%mu=4e-7*pi;
mu=1;
gam=26.7519E+7;
M=zeros(length(t1), length(t2));
for i=1:length(t1)
	M(i,:)=sin(w*t1(i)+phi)*exp(complex(0,w*t2+theta)).*exp(-D*Mo*t2*complex(0,1)*cos(w*t1(i)+phi));
end
figure(3);
%plot(real(M));
M=abs(fftshift(fft2(M)));
%M=real(M);
%pcolor(M); shading interp;
%plot(sum(M,2))
plot(M)
plot(imag(exp(-D*Mo*t2*complex(0,1)*cos(w*tau+phi))))

figure(4)
lim=100; rr=-lim:lim;
M=zeros(size(t2));
for i=-lim:lim
	M=M+(complex(0,-1)^i)*besselj(i,D*Mo*t2)*exp(complex(0,-1)*i*(w*tau+phi)).*exp(complex(0,1)*(w*t2));
end
plot(real(M))
