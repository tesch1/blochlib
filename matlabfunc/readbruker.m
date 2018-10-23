
function [dat, varargout]=readbruker(dirr)
% [data, {sweep, is2D}]=readbruker('dir')
%
% author:: Bo Blanton UC Berkeley, Dept of Chemistry
% email:: magneto@dirac.cchem.berkeley.edu
% more info:: http://waugh.cchem.berkeley.edu/blochlib/
% last modified:: 10.20.02
%
%   reads a bruker file directory to get the FID(s)
%   it requires the normal XWINNMR parameters file
%   'acqu' and the 'fid' if 1D and 'ser' if 2D
%   returns the data in 'data' and the sweep width in 'sweep'


is2D=0;

dat=[];

acquF=[dirr '/acqu'];
tm=fopen(acquF);
if(tm == -1)
	error(['Directory, ' dirr 'does not contain an "acqu" file']);
	return;
else
	fclose(tm);
end

acquF2D=[dirr '/acqu2'];
tm=fopen(acquF2D);
if(tm ~= -1)
	is2D=1;
else
	fclose(tm);
end

td	=	brxget(acquF, 'TD', 'n');
sw	=	brxget(acquF, 'SW_h', 'n');
decim =128;

tacq=td/(2*sw);		%total acquisition time

td2=1;
if(is2D==1)
	td2	=	brxget(acquF2D, 'TD', 'n');
end

%looking at the particular type of data

if (is2D==0)
	datafile=[dirr '/fid']; 	%is2D=0 for fid
else
	datafile=[dirr '/ser'];		%is2D=1 for ser
end

%test to see that the datafile exsists
FIDDATA=fopen(datafile,'r','b');
if(FIDDATA==-1)
	error(['Directory, ' dirr 'does not contain an "fid" or "ser" file']);
	return;
end	

%the fid sizes
size1D=td/2;
size2D=td2;

dat=zeros(size1D, size2D);

%
for cnt = 1:size2D;
	
	%'int32' is the identifier of the type of data
	
	%brfid contains the td points of the fid.  It is a row vector.

	brfid=fread(FIDDATA,td,'int32');

	%'reshape' sorts data into real and imaginary parts 
	%reshape can take every other point, first into one rove
	%next into another row of a 2 by td/2 matrix

 	%RESHAPE change size.
	%    RESHAPE(X,M,N) returns the M-by-N
	%matrix whose elements
	%    are taken columnwise from X.  An
	%error results if X does
	%    not have M*N elements.

	brfid=reshape(brfid,2,length(brfid)/2);

	%now the two rows of the matrix brfid2 are
	%the real and imag parts of the fid.  The
	%parts can be individually accessed using
	%the colon ':'

 	%realfid=brfid(1,:);

 	%imagfid=brfid(2,:);

	%formation of the complex vector

	cplxfid=brfid(1,:)+i*brfid(2,:);

	%cplxfid has a form real+i*imag

	%the correction for dc offset
	%the first 75 points of the FID in doing the DC offset
	%correction are ignored.
	%these points are artificial if the digital filter is on.

	cplxfid=cplxfid-mean(cplxfid(75:end));

	%the data is placed in a column vector by taking
	%the transpose with the .' command.

	cplxfid = cplxfid.';


	%Bruker data acquired with the digital filter on has
	%around 70 'bad points' at the beginning of the data

	%These points will be removed from there. 
	%normally they should be moved to the end of the fid
	%but in case of our calculation it is fine to take them 
	%out

	%The number of points to be removed are calculated from 
	%decim parameter in acqus file:	%
	% points = (70.5 - 15.5/DECIM)  for DECIM a power of 2
	% points = (185/3 - 15.5/DECIM)  for DECIM not a power of 2
	
	%determine if decim is a power of 2

 
	test=log2(decim);
	test=test-fix(test);

	%test=0 if decim is a power of 2
	%calculate how many points to remove

	if (test==0)
		points=(70.5-15.5/decim);
	else
		points=(185/3-15.5/decim);
	end

	nskip=fix(points);

	cplxfid=cplxfid(nskip:end);

	%ncplxfid=cplxfid(1:points); 	%removed points

	%it is stored as the cnt-th column of a matrix with


	dat(:,cnt)=cplxfid;

	%nfidmat(:,cnt)=ncplxfid; 	%removed points

end

if(nargout>=2) varargout(1) = {sw}; end
if(nargout>=3) varargout(2) = {is2D}; end

return;

%------------------------------------------
%sub function that reads the acq file
function value = brxget(inname, srcstring, type)
% function value = brxget(inname,srcstring, type)
%
% reads XWINNMR acquisition file inname and returns 
% the value of the parameter srcstring
% 
% input:
%   inname--> is a string with the filename of the file
%   srcstring--> is a string that is the name of 
%      the parameter to be found
%   type--> is either 'n' or 't' specifying if we are
%      searching for a number or text
% output:
%   value --> is the value of the parameter
%     if the parameter is NOT found 
%     value=pi (just becuase the file will not likely
%               have 'pi' as a parameter somewhere)


	eval(['h = fopen(''' inname ''');'])
	value=pi;
	noteof = 1; % switch when we reach end of file
	a = 0;
	findme = [srcstring, '= '];

	while noteof % read in the file line by line    
		line = fgetl(h); % read in next line of the file
		if ischar(line)
			a = findstr(line, findme);
			if isempty(a)~=1
				if type == 'n'
					value = str2num(line(a+length(findme):end));
				else
					value = line(a+length(findme):end); 
				end
			end
		else
			noteof = 0;
		end 
	end
	

return;

%------------------------------------------
