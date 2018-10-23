function [dat, varargout]=readspinsight(dirr)
% [data, {sweep, is2D}]=readspinsight('dir')
%
% author:: Bo Blanton UC Berkeley, Dept of Chemistry
% email:: magneto@dirac.cchem.berkeley.edu
% more info:: http://waugh.cchem.berkeley.edu/blochlib/
% last modified:: 10.20.02
%
%   reads a Spin Sight file directory to get the FID(s)
%   it requires the normal Spin Sight parameters file
%   'acq' and the 'data' for both 1D and 2D
%
%   returns the data in 'data' 
%
% optional ouputs are sweep width in 'sweep'
% and 'is2D' will be '1' if the data is 2D
% and 0 if the data is 1D

is2D=0;

dat=[];

acquF=[dirr '/acq'];
tm=fopen(acquF);
if(tm == -1)
	error(['Directory, ' dirr 'does not contain an "acq" file']);
	return;
else
	fclose(tm);
end

%get the valid parameters in the files
td =	spinsightget(acquF, 'al', 'n');
dw =	spinsightget(acquF, 'dw', 'n'); %dwell time
sw=1.0/dw; 
	
td2=spinsightget(acquF, 'al2', 'n');
if(td2==pi)
	is2D=0;
	td2=1;
else
	is2D=1;
end

%looking at the particular type of data

datafile=[dirr '/data']; 	

%test to see that the datafile exsists
FIDDATA=fopen(datafile,'r','b');
if(FIDDATA==-1)
	error(['Directory, "' dirr '" does not contain an "data" file']);
	return;
end	

%the fid sizes
size1D=td;
size2D=td2;

dat=zeros(size2D, size1D);
a=[];
count=0;
if size1D~=pi | size2D~=1
	[a, count] = fread(FIDDATA, 2*size1D*size2D, 'int32');
else
	fseek(FIDDATA, 0,1);
	len=ftell(FIDDATA);
	fseek(FIDDATA, 0, -1);
	size1D=len/4/2;
	[a, count] = fread(FIDDATA, 2*size1D*size2D, 'int32');
end
disp(ferror(FIDDATA));
fclose(FIDDATA);
for tel=1:size2D
    dat(tel, 1:(size1D)) = (a( ((tel-1)*size1D + 1) : ((tel)*size1D)) + sqrt(-1)*a( ((tel-1)*size1D + size2D*size1D + 1) : ((tel)*size1D + size2D*size1D)) )';

%	dat(tel, 1:(size2D)) = (a( ((tel-1)*size2D + 1) : ((tel)*size2D)) + sqrt(-1)*a( ((tel-1)*size2D + size1D*size2D + 1) : ((tel)*size2D + size1D*size2D)) )';
end;
if(nargout>=2) varargout(1) = {sw}; end
if(nargout>=3) varargout(2) = {is2D}; end
return;



%------------------------------------------
%sub function that reads the acq file
function value = spinsightget(inname, srcstring, type)
% function value = spinsightget(inname,srcstring, type)
%
% reads spin sight acquisition file inname and returns 
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

	noteof = 1; % switch when we reach end of file
	a = 0;
	findme = [srcstring, '='];
	value=pi;	
	while noteof % read in the file line by line    
		line = fgetl(h); % read in next line of the file
		if ischar(line)
			a = findstr(line, findme);
			if isempty(a)~=1
				if type == 'n'
					value = str2double(line(a+length(findme):length(line)));
					if isnan(value)
						value = str2double(line(a+length(findme):length(line)-1));
					end
				else
					value = line(a+length(findme):length(line)); 
				end
			end
		else
			noteof = 0;
		end 
	end

return;

%------------------------------------------
