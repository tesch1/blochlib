

MAXITER=10000;

R=[2 5 20 24 31 32 36 48 64 73 96 97 127 128 129 163 191 192 229 255 256 257 319 320 321 417 479 480 511 512 ];
matlab_prop=1:length(R);

mytime=cputime;
aftertime=0;
holdt=1;
for n=1:length(R);
	A=complex(rand(R(n),R(n)), rand(R(n),R(n)));
	B=complex(rand(R(n),R(n)), rand(R(n),R(n)));
	for i=1:MAXITER
		nowt=cputime;
		C=A*B*A';
		aftertime=aftertime+(cputime-nowt);
		if (aftertime*1e6)>2.0
			holdt=i;
			break;		
		end
		holdt=MAXITER;
	end
	totaltim=(aftertime)/holdt; 
	
% each mat mul is M^3 (and there are 2)
% a complex mat mul is the same a 4 non-complex mat muls
% the adjint op is a M operations
% finally a double is 2 floats (for a MFLOPS count)
	Nflops=( 2.0 *4.0*2.0 * R(n) * R(n) * R(n) +R(n))/ 1.0e6;
	sprintf('For %d x %d complex cfoo*cfoo*adjoint(cfoo): %g micro-sec: MFLOPS %g', R(n),R(n),totaltim*1e6, Nflops/totaltim)
	matlab_prop(n)=Nflops/totaltim;
end

matlab_mul=1:length(R);

mytime=cputime;
aftertime=0;
holdt=1;
for n=1:length(R);
	A=complex(rand(R(n),R(n)), rand(R(n),R(n)));
	B=complex(rand(R(n),R(n)), rand(R(n),R(n)));
	for i=1:MAXITER
		nowt=cputime;
		C=A*B;
		aftertime=aftertime+(cputime-nowt);
		if (aftertime*1e6)>2.0
			holdt=i;
			break;		
		end
		holdt=MAXITER;
	end
	totaltim=(aftertime)/holdt; 
	
% each mat mul is M^3 
% a complex mat mul is the same a 4 non-complex mat muls
% finally a double is 2 floats (for a MFLOPS count)
	Nflops=( 4.0*2.0 * R(n) * R(n) * R(n))/ 1.0e6;
	sprintf('For %d x %d complex cfoo*cfoo: %g micro-sec: MFLOPS %g', R(n),R(n),totaltim*1e6, Nflops/totaltim)
	matlab_mul(n)=Nflops/totaltim;
end

