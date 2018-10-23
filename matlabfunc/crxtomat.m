function [is1D, datout]=crxtomat(dat)

% [is1D, datout]=crxtomat(dat)
%
%takes a chemmagnetics ASCII data file and converts
%it to a matlab array...
%
% author:: Bo Blanton UC Berkeley, Dept of Chemistry
% email:: magneto@dirac.cchem.berkeley.edu
% more info:: http://waugh.cchem.berkeley.edu/blochlib/
% last modified:: 10.20.02
%
% 'dat' the previously loaded ascii data
% with the following format
% 
% point(dim1) point(dim2) fidamp
%     0            0         # 
%     1            0         #
%    ...          ...       ... 
%     0            0         #
%     1            0         #
%
%    the first set of numbers is the real parts
%    the second set is for the imag parts
%
%  1) load the ascii data file first
%  2) send the variable to this function.
%
%    >>load mydat.asc;
%    >>matdat=crxtomat(mydat);
%    >>[is1D matdat]=crxtomat(mydat);

%first we need to find the 'length' of the first dim
%and the length of the second dim

if size(dat, 2) == 3
	maxsize=size(dat,1);
	maxx=dat(maxsize, 1)+1;
	maxy=dat(maxsize, 2)+1;
	mx=maxx*maxy;
	tr=dat(1:(mx),3);
	ti=dat((mx+1):(2*mx),3);
	trm=reshape(tr, maxx, maxy);
	tim=reshape(ti, maxx, maxy);

	datout=trm+complex(0,1)*tim;
	is1D=0;
	return
elseif size(dat, 2) == 2
	maxsize=size(dat,1);
	maxx=dat(maxsize, 1)+1;
	mx=maxx;
	tr=dat(1:(mx),2);
	ti=dat((mx+1):(2*mx),2);

	datout=tr+complex(0,1)*ti;
	is1D=1;
	return
end	

   
