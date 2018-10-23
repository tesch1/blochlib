function dat = phaser(data, center,ph0,ph1);
% dat=phaser(data,center, ph0, ph1)
%
% author:: Bo Blanton UC Berkeley, Dept of Chemistry
% email:: magneto@dirac.cchem.berkeley.edu
% more info:: http://waugh.cchem.berkeley.edu/blochlib/
% last modified:: 10.20.02
%
% A phasing function that will work on both 1D and 2D data sets
% data--> the input data
% center-->the 'center' point to perform the phasing around
% ph0 --> zero order phase correction
% ph1--> first order phase corection
%

if nargin < 3;
  ph1 = 0;
end;
if nargin < 2;
  ph0 = 0;
end;
[m n]=size(data);
dat=zeros(m,n);
center=abs(center);
i=sqrt(-1);

if m==1
	z=length(data(1,:));
	if center>z
		center=z;
	end
	z=(-(center)):1:(z-center-1);
	z=z/(length(z)-1);
	dat=data.* exp(i*(ph0/180*pi+ph1/180*pi*z));

elseif n==1
	z=length(data(:,1));
	if center>z
		center=z;
	end
	z=(-(center)):1:(z-center-1);
	z=z/(length(z)-1);
	dat=data.*exp(i*(ph0/180*pi+ph1/180*pi*z)).';
else
	z=length(data(1,:));
	if center>z
		center=z;
	end
	z=(-(center)):1:(z-center-1);
	z=z/(length(z)-1);
	for r=1:m;
		dat(r,:)=data(r,:).* exp(i*(ph0/180*pi+ph1/180*pi*z));
	end;
end




