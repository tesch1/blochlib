function [x,y,z,v]=makemesh(data,xd,yd,zd, slices)

% [x,y,z,v]=makemesh(data, Nx, Ny, Nz, [slices])
%
%a function that takes a x,y,z, SCALAR data set and makes 
% nice 3D gaphs from the data...it also returns the 'mesh grids'
% needed to be generated in order to create the graphs
%
% data -> the input data in the following format
%          x1 y1 z1 data
%          x2 y1 z1 data
%          x3 y1 z1 data
%            .....
%          xNx y1 z1 data
%          x1 y2 z1 data
% Nx-> the number of point of the X dimension (this is NOT the file length)
% Ny-> number of points along the Y dimension (this is NOT the file length)
% Nz-> number of points along the Z dimension (this in NOT the file length)
%
%  the matrix sizes of [x,y,z,v] will be Nx x Ny x Nz

ct=1;
if size(data, 2)==6
	ct=3; 
end;
x=squeeze(data(1:xd,1));
y=reshape(data(:,2), xd,yd,zd);
y=squeeze(y(1,:,1));
z=reshape(data(:,3), xd,yd,zd);
z=squeeze(z(1,1,:));
v={};
if nargin<5
	slices=1;
end
v=reshape(data(:,3+slices), xd,yd,zd);
%for curPt=1:ct
	figure(slices);

	%subplot(1,2,1);
	hold on;
	maxx=max(x(:));
	maxy=max(y(:));
	maxz=max(z(:));


	minx=min(x(:));
	miny=min(y(:));
	minz=min(z(:));
	[xx yy zz]=meshgrid(x,y,z);

	hslice = surf(linspace(minx,maxx,100),...
	    linspace(minz,maxz,100),zeros(100));

	rotate(hslice,[-1,0,0],-70)
	xd = get(hslice,'XData');
	yd = get(hslice,'YData');
	zd = get(hslice,'ZData');
	delete(hslice)
	h = slice(x,y,z,v,xd,yd,zd);
	set(h,'FaceColor','interp',...
	    'EdgeColor','none',...
	    'DiffuseStrength',.8)




	hx = slice(xx,yy,zz,v,maxx,[],[]);
	set(hx,'FaceColor','interp','EdgeColor','none')

	hy = slice(xx,yy,zz,v,[],maxy,[]);
	set(hy,'FaceColor','interp','EdgeColor','none')

	hz = slice(xx,yy,zz,v,[],[],minz);
	set(hz,'FaceColor','interp','EdgeColor','none')

	%daspect([1,1,1.5])
	axis tight
	box on
	view(-38.5,16)
	camzoom(1)
	%camproj perspective
	lightangle(-45,45)

	set(gcf,'Renderer','zbuffer')
	hold off;
	%subplot(1,2,2);
	%plot3(data(:,1), data(:,2), data(:,3), 'k.');		
%end

