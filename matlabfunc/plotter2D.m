function h=plotter2D(FID, num, NC,sw1, sw2)
% A plotter tool for 2D data...does no processing, just plots the 
% data in various plot types
% use 'ffter' to perform signal processing to the data
%
% author:: Bo Blanton UC Berkeley, Dept of Chemistry
% email:: magneto@dirac.cchem.berkeley.edu
% more info:: http://waugh.cchem.berkeley.edu/blochlib/
% last modified:: 10.20.02
%
%plotter2D(FID,num, {NC, sw1, sw2})
% FID--> your 2D data set
% num--> the plot id
% 1--3 ARE VALID ONLY FOR 1D DATA SETS ONLY
%  1=simple 1D plot (real)
%  2=simple 1D plot (imag)
%  3=simple 1D plot (abs)
% 11--83 ARE VALID FOR 2D DATA SETS ONLY
%  11=pcolor (real)
%  12=pcolor (imag)
%  13=pcolor (abs)
%  21=the collaped sum of the 1st dimendsion (real)
%  22=the collaped sum of the 1st dimendsion (imag)
%  23=the collaped sum of the 1st dimendsion (abs)
%  31=the collaped sum of the 2nd dimendsion (real)
%  32=the collaped sum of the 2nd dimendsion (imag)
%  33=the collaped sum of the 2nd dimendsion (abs)
%  41=the FID as a surface (real)
%  42=the FID as a surface (imag)
%  43=the FID as a surface (abs)
%  51=20 level contour (real)
%  52=20 level contour (imag)
%  53=20 level contour (abs)
%  61=a stacked plot (real part) (like the 2D display in the CMX software)
%  62=a stacked plot (imag part) (like the 2D display in the CMX software)
%  63=a stacked plot (abs) (like the 2D display in the CMX software)
%  71=a concatonated (real) 1D plot of each slice side by side
%  72=a concatonated (imag) 1D plot of each slice side by side
%  73=a concatonated (abs) 1D plot of each slice side by side
%  81=plot3 (real) a stacked splot in 3D
%  82=plot3 (imag) a stacked splot in 3D
%  83=plot3 (abs) a stacked splot in 3D
% NC --> valid for contour plots (optional, default is 20)
% sw1--> sweepwidth for the first dim (optional)
% sw2--> sweepwidth for the second dim (optional)

if(nargin==1 | nargin==0)
	error(' must specify FID and a plot number' )
end

NConts=20;
if nargin>=3
	NConts=NC;
end

if num<10
	xsiz=size(FID, 1);
	ysiz=size(FID,2);
	if ysiz>xsiz & (xsiz==1 | ysiz==1)
		FID=FID.';
	end
end

xlab='Second Dim';
ylab='First Dim';
DASpts=size(FID,2);
MASpts=size(FID,1);
X=1:MASpts;
Y=1:DASpts;


if(nargin>=4)
	if sw1 ~= 0  & MASpts>1
		if length(sw1) ~= MASpts
			if length(sw1)>1 
				MASsweep=max(sw1)*2;
			else
				MASsweep=sw1;
			end
			X=(-MASsweep/2:MASsweep/(MASpts-1):MASsweep/2);
		else
			X=sw1;
		end
	end	
end

if(nargin>=5)
	if sw2~=0 & DASpts>1
		if length(sw2) ~= DASpts
			if length(sw2)>1
				wr=max(sw2)*2;
			else
				wr=sw2;
			end
			Y=(-wr/2:wr/(DASpts-1):wr/2);
		else
			Y=sw2;
		end
	end
end
minshowMAS=min(Y);
maxshowMAS=max(Y);


if num>10
	X=fliplr(X);

	if num==21 | num==31
		DASdat=real(sum(FID,1));
		MASdat=real(sum(FID,2));
		Y=fliplr(Y);
	elseif num==22 | num==32
		DASdat=imag(sum(FID,1));
		MASdat=imag(sum(FID,2));
		Y=fliplr(Y);
	elseif num==23 | num==33
		DASdat=abs(sum(FID,1));
		MASdat=abs(sum(FID,2));
		Y=fliplr(Y);
	end
	
	lM=1;
	rM=MASpts;
	lD=1;
	rD=DASpts;

	for i=1:MASpts
		if(X(i)<=minshowMAS)
			lM=i;
		elseif(X(i)>=maxshowMAS)
			rM=i;
			break;
		end
	end
end

if num==1
	plot(X, real(FID));
elseif num==2
	plot(X, imag(FID));
elseif num==3
	plot(X, abs(FID));
elseif num==81
	[XX,YY]=meshgrid(Y,X);
	plot3(XX, YY, real(FID),'k');
elseif num==82	
	[XX,YY]=meshgrid(Y,X);
	plot3(XX, YY, imag(FID),'k');
elseif num==83   
	[XX,YY]=meshgrid(Y,X);
	plot3(XX, YY, abs(FID),'k');
elseif num==11	
	[XX,YY]=meshgrid(Y,X);
	pcolor(XX,YY,real(FID));shading 'interp'
elseif num==12	
	[XX,YY]=meshgrid(Y,X);
	pcolor(XX,YY,imag(FID));shading 'interp'
elseif num==13   
	[XX,YY]=meshgrid(Y,X);
	pcolor(XX,YY,abs(FID)); shading 'interp'
elseif num==21 | num==22 | num==23  
	plot(X,MASdat);
elseif num==31 | num==32 | num==33 
	plot(Y,DASdat);
elseif num==41	
	pp=real(FID);  
	j=surf(Y,X,pp);   
	set(j,'EdgeColor','none');   
	set(gca, 'ZTick', []);	
	xlabel(xlab);	
	ylabel(ylab);   
	brighten(.7);   
	lighting phong;  
	for i =1:DASpts		
		maxxs(i)=real(max(pp(:,i)));  
	end   
	mightymax=max(maxxs);   
	material metal;
	axis([-Inf Inf -Inf Inf -Inf Inf]);  
elseif num==42	
	pp=imag(FID);  
	j=surf(Y,X,pp);   
	set(j,'EdgeColor','none');   
	set(gca, 'ZTick', []);	
	xlabel(xlab);	
	ylabel(ylab);   
	brighten(.7);   
	lighting phong;  
	for i =1:DASpts		
		maxxs(i)=real(max(pp(:,i)));  
	end   
	mightymax=max(maxxs);   
	material metal;
	axis([-Inf Inf -Inf Inf -Inf Inf]);  
elseif num==43	
	pp=abs(FID);  
	j=surf(Y,X,pp);   
	set(j,'EdgeColor','none');   
	set(gca, 'ZTick', []);	
	xlabel(xlab);	
	ylabel(ylab);   
	brighten(.7);   
	lighting phong;  
	for i =1:DASpts		
		maxxs(i)=real(max(pp(:,i)));  
	end   
	mightymax=max(maxxs);   
	material metal;
	axis([-Inf Inf -Inf Inf -Inf Inf]);  
elseif num==51
	[xx,yy]=meshgrid(Y,X);
	contour(xx,yy,real(FID),NConts);
elseif num==52
	[xx,yy]=meshgrid(Y,X);
	contour(xx,yy,imag(FID),NConts);
elseif num==53
	[xx,yy]=meshgrid(Y,X);
	contour(xx,yy,abs(FID),NConts);

elseif num==61
	zz=real(FID);	
	surf(Y,X,zz, 'MeshStyle', 'column', 'FaceColor', 'white','Clipping', 'on' );
	set(gca, 'ZTick', []);	
	axis([minshowMAS maxshowMAS -Inf Inf -Inf Inf]);
	view(90,78);
	
elseif num==62
	zz=imag(FID);	
	surf(Y,X,zz, 'MeshStyle', 'column', 'FaceColor', 'white','Clipping', 'on' );
	set(gca, 'ZTick', []);	
	axis([minshowMAS maxshowMAS -Inf Inf -Inf Inf]);
	view(90,78);
elseif num==63
	zz=abs(FID);	
	surf(Y,X,zz, 'MeshStyle', 'column', 'FaceColor', 'white','Clipping', 'on' );
	set(gca, 'ZTick', []);	
	xlabel(xlab);	
	ylabel(ylab);   
	axis([minshowMAS maxshowMAS -Inf Inf -Inf Inf]);
	view(90,78);
	
elseif num==71
	zz=reshape(FID, DASpts*MASpts, 1);
	plot(real(zz));
elseif num==72
	zz=reshape(FID, DASpts*MASpts, 1);
	plot(imag(zz));

elseif num==73
	zz=reshape(FID, DASpts*MASpts, 1);
	plot(abs(zz));
	
else	
	disp('not a valid number');
end	
