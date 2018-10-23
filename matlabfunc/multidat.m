
function dat=multidat(fname, initf,endf)

%multiplot(fname, initf, endf)
%
% author:: Bo Blanton UC Berkeley, Dept of Chemistry
% email:: magneto@dirac.cchem.berkeley.edu
% more info:: http://waugh.cchem.berkeley.edu/blochlib/
% last modified:: 10.20.02
%
% takes a base file name for a data series in separate files
% and concatonate them into one 2-D data set
% then uses 'nmrplot' to plot the result
% 
% fname--> the base file name
%  files should be in a series like "data0", "data1", "data2"...
%  it is ASSUMED that all the files have the same data lengths
% files--> the number of files...
% Returns the matrix of data
%

datlen=0;
dat=0;
ct=1;
for i=initf:endf
	loo=strcat(fname,sprintf('%d', i));
	load(loo);
	sname=eval(loo);
	if i==initf
		datlen=size(sname,1);
		dat=zeros(datlen, abs(endf-initf));
	end
	dat(:,ct)=complex(sname(:,5), sname(:,6));
	ct=ct+1;
end

return 
	



