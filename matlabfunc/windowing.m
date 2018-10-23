function fid=windowing(fidin, tty, p1,p2,p3)

%fid=windowing(fidin, type, {p1, p1, p3})
%
% author:: Bo Blanton UC Berkeley, Dept of Chemistry
% email:: magneto@dirac.cchem.berkeley.edu
% more info:: http://waugh.cchem.berkeley.edu/blochlib/
% last modified:: 10.20.02
%
%A Basic windowing function with a few 'types' of windowing
% type=0-->Do NOTHING
% type=1--> Exponential (requires 'p1' as a parameters) 
%           does exp(p1)*fidin
% type=2--> Gaussain (requires 'p1' as the 'offset' and 'p2' as the sigma) 
%           does exp((i-p1)^2/2*p2*p2)*fidin
% type=3--> Hanning (requires p1 as an 'offset')
%           does (0.5+0.5*cos(Pi/size(fidin)*(i-p1)))*fidin
% type=4--> Hamming (requires 'p1' as an 'offset')
%           does (0.54 + 0.46*cos(Pi/size(fidin)*(i-p1)))*fidin


if nargin<2
	disp('ERROR:: Windowing must contain at least 2 arguments');
else
	st1=size(fidin,1);
	st2=size(fidin,2);
	if st1~=1 & st2~=1
		disp('ERROR::Windowing works only for VECTORs, not matrices');
		fid=fidin;
		return;
	end
	len=length(fidin);
	fid=zeros(size(fidin));
	if tty==0
		fid=fidin;
		return;
	elseif tty==1
		if nargin<3
			disp('ERROR::Exponential windows must have 3 arguments');
			fid=fidin;
			return 
		else
			for i=0:(len-1)
				fid(i+1)=exp(-abs(p1)*i/len).*fidin(i+1);
			end
			return;
		end
	elseif tty==2
		if nargin<4
			disp('ERROR::Gaussian windows must have 4 arguments');
			fid=fidin;
			return;
		else
			denom = 2*p2*p2;
			for i=1:len			    
			    x = (i/len-p1);			
			    fid(i) = exp(-x*x/denom)*fidin(i);
			end
			return;
		end
	elseif tty==3
		if nargin<3
			disp('ERROR::Hanning windows must have 3 arguments');
			fid=fidin;
			return
		else
			fact=pi/len;
			for i=1:len
				x=(i-p1);
				fid(i)=(0.5+0.5*cos(fact*x))*fidin(i);
			end
			return;
		end
	elseif tty==4
		if nargin<3
			disp('ERROR::Hamming windows must have 3 arguments');
			fid=fidin;
			return;
		else
			fact=pi/len;
			for i=1:len
				x=(i-p1);
				fid(i)=(0.54+0.46*cos(fact*x))*fidin(i);
			end
			return;
		end
	else
		disp('ERROR:: windowing option not available');
		fid=fidin;
		return ;
	end
end

		

