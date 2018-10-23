

MAXITER=100;

n=36;
cfoo=complex(rand(n,n), rand(n,n));
tmpfoo=cfoo;
rfoo=rand(n,n);
assfoo=[];

mytime=cputime;
aftertime=0;
for i=1:MAXITER
	nowt=cputime;
	assfoo=rfoo*rfoo;
	aftertime=aftertime+(cputime-nowt);
	i=i+1;
end
totaltim=(aftertime)/MAXITER; 
sprintf('For %d x %d real cfoo*cfoo: %g micro-sec: MFLOPS %g', n,n,totaltim*1e6, 2.0 * n * n * n / totaltim/1e6)

mytime=cputime;
aftertime=0;
for i=1:MAXITER
	nowt=cputime;
	assfoo=cfoo*cfoo;
	aftertime=aftertime+(cputime-nowt);
	i=i+1;
end
totaltim=(aftertime)/MAXITER; 
sprintf('For %d x %d complex cfoo*cfoo: %g micro-sec: MFLOPS %g', n,n,totaltim*1e6, 4.0 *  n * n * n / totaltim/1e6)

mytime=cputime;
for i=1:MAXITER
	assfoo=cfoo*cfoo*cfoo';
	i=i+1;
end
aftertime=cputime;
totaltim=(aftertime-mytime)/MAXITER; 
sprintf('For %d x %d complex cfoo"*cfoo*cfoo: %g micro-sec: MFLOPS %g', n,n,totaltim*1e6,4.0 *  n * n * n*3 / totaltim/1e6)

mytime=cputime;
for i=1:MAXITER
	assfoo=expm(cfoo);
	i=i+1;
end
aftertime=cputime;
totaltim=(aftertime-mytime)/MAXITER; 
 
sprintf('For %d x %d complex matrix diagonalization EXPM: %g micro-sec', n,n,totaltim*1e6)


mytime=cputime;
for i=1:MAXITER
	assfoo=expm1(cfoo);
	i=i+1;
end
aftertime=cputime;
totaltim=(aftertime-mytime)/MAXITER; 
 
sprintf('For %d x %d complex matrix diagonalization EXPM1: %g micro-sec', n,n,totaltim*1e6)

mytime=cputime;
for i=1:MAXITER
	assfoo=expm2(cfoo);
	i=i+1;
end
aftertime=cputime;
totaltim=(aftertime-mytime)/MAXITER; 
 
sprintf('For %d x %d complex matrix diagonalization EXPM2: %g micro-sec', n,n,totaltim*1e6)

mytime=cputime;
for i=1:MAXITER
	assfoo=expm3(cfoo);
	i=i+1;
end
aftertime=cputime;
totaltim=(aftertime-mytime)/MAXITER; 
 
sprintf('For %d x %d complex matrix diagonalization EXPM3: %g micro-sec', n,n,totaltim*1e6)


mytime=cputime;
for i=1:MAXITER
	assfoo=expm(rfoo);
	i=i+1;
end
aftertime=cputime;
totaltim=(aftertime-mytime)/MAXITER; 
 
sprintf('For %d x %d real matrix diagonalization: %g micro-sec', n,n,totaltim*1e6)
