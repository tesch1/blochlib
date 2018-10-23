function ffts=ffter(fids, appo1, appo2, dcoff,zfil,fftyn, windopt, windval1, windval2)

%ffts=ffter(fids, appo1, appo2, dcoff,zfil,fftyn,{windopt, windval1, windval2})
%
%little function that preforms a
% dcoffset correction and fourier transfrom
% to a matrix of fids
% 
% author:: Bo Blanton UC Berkeley, Dept of Chemistry
% email:: magneto@dirac.cchem.berkeley.edu
% more info:: http://waugh.cchem.berkeley.edu/blochlib/
% last modified:: 10.20.02
%
% appo1--> appodization along the t1 axis
% appo2--> appodization along the t2 axis
% dcoff--> correct dc offsets
% zfill--> zero fill the data
% fftty--> 
%   if 1 scalar value
%	0=NO fft, 
%	1=full FFT, 
%	2=FFT real part Only
%	3=FFT imag part only
%  if a 2x1 array
%	The same options as above except the first element 
% 	is for the first dimension
%	the second element for the second dimension
% {windopt} --> extra option for windows type (see 'windowing.m')
% {windval1} --> extra option for windowing
% {windval2) --> extra option for windowing


st1=size(fids,1);
st2=size(fids,2);

if st2==1
	fids=fids';
	st1=size(fids,1);
	st2=size(fids,2);
end
	
%ffts=zeros(st1*2,st2*2);
if(zfil==1)
	if st1==1
		newst1=st1;
		newst2=2*st2;
	else
		newst1=2*st1;
		newst2=2*st2;
	end
		
else
	newst1=st1;
	newst2=st2;
end

appd1=exp(-abs(appo2)*(0:(newst1-1))/newst1);
appd2=exp(-abs(appo1)*(0:(st2-1))/st2);

fftynd1=0;
fftynd2=0;
if max(size(fftyn))>=2
	fftynd1=fftyn(1);
	fftynd2=fftyn(2);
else
	fftynd1=fftyn;
	fftynd2=0;
end

ffts=zeros(newst1,newst2);

if st1==1
	if(dcoff==1)
		loo=mean(fids(ceil(st2/2):ceil(st2)));
	else
		loo=0;
	end
	

	fids=fids-loo;
	fids=fids.*appd2;  
	if nargin>8
		fids=windowing(fids, windopt, windval1, windval2);
	end

	fidd=zeros(1,newst2);
	fidd(1:st2)=fids(1:st2);   
	if(fftynd1==1 )
		ffts=fftshift(fft(fidd));
	elseif(fftynd1==2)
		ffts=fftshift(fft(real(fidd)));
	elseif(fftynd1==3)
		ffts=fftshift(fft(imag(fidd)));
	else
		ffts=fidd;
	end
else
	for i=1:st1
		if(dcoff==1)
			loo=mean(fids(i,ceil(st2/2):ceil(st2)));
		else
			loo=0;
		end
		

		fids(i,:)=fids(i,:)-loo;
		fids(i,:)=fids(i,:).*appd2;  
		if nargin>8
			fids(i)=windowing(fids(i), windopt, windval1, windval2);
		end
		fidd=zeros(newst2,1);
		fidd(1:st2)=fids(i,1:st2);   
		
		if(fftynd1==1)
			ffts(i,1:newst2)=fftshift(fft(fidd)).';
		elseif(fftynd1==2)
			ffts(i,1:newst2)=fftshift(fft(real(fidd))).';
		elseif(fftynd1==3)
			ffts(i,1:newst2)=fftshift(fft(imag(fidd))).';
		else
			ffts(i,1:newst2)=fidd.';
		end
	end
end
st2=st2*2;

if st1>1
	for i=1:newst2     
		
   		if(dcoff==1)
      			loo=mean(ffts(ceil(st1/2):ceil(st1),i));
      		else
      			loo=0;
      		end
  		
  		ffts(:,i)=ffts(:,i).*appd1.';
		if nargin>8
			fids(i)=windowing(fids(i), windopt, windval1, windval2);
		end

      		if(fftynd2==1)
      			ffts(:,i)=fftshift(fft(ffts(:,i)-loo));
      		elseif(fftynd2==2)
      			ffts(:,i)=fftshift(fft(real(ffts(:,i)-loo)));
		elseif(fftynd2==3)
      			ffts(:,i)=fftshift(fft(imag(ffts(:,i)-loo)));
		else
      			ffts(:,i)=ffts(:,i)-loo;
      		end
      		
   	end
end
return
