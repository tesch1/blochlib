
function [XX,YY,ZZ,B0,B1,B2]=plotmag(indata, which, isoslices, numlines, view, plotcoil)

% [X,Y,Z, Bx,By,Bz]=plotmag(data, which, [ isoslices, numlines,view, plotcoil])
%
% author:: Bo Blanton UC Berkeley, Dept of Chemistry
% email:: magneto@dirac.cchem.berkeley.edu
% more info:: http://waugh.cchem.berkeley.edu/blochlib/
% last modified:: 10.20.02
%
%  Plots various types of plots for Magnetic Field data as generated by 
%  'biot.h' and 'biot.cc' in the 'Blochlib' pakage 
%  (http://waugh.cchem.berkeley.edu/blochlib/)
%
% 'data'--> is the data file to be read in
% 'which' --> is the plot type
%        =1 --> a simple plot3 of the coil and the grid points overlaid
%        =2 --> a 3x3 bunch of contour plots that allow you to move through them
%        =3 --> a 3x1 3D rendered bunch of iossurfaces
%        =4 --> a 3x1 3D contour and pcolor surfaces (uses 'alpha' if matlab 6)
% 'isoslices' --> is only valid for option which=3..the number of silces to include in the 3D renders
% 'numlines'--> only for options which=3 and 4...determins divisions the stream lines
% 'view' --> 'x', 'y', 'z', or 'a' (a==all)
% 'plotcoil' --> 1 or 0...add the coil to the plots (for which==4 ONLY)
% it also returns the 3D grid strucutures for the grid and the Magnetic field
%

load(indata);

vr=version;
B0=reshape(B0,N0,N1,N2);
B1=reshape(B1,N0,N1,N2);
B2=reshape(B2,N0,N1,N2);

vx=grid0(1:N0);
vy=grid1(1:N0:N1*N0);
ZZ=reshape(grid2, N0,N1,N2);
vz=ZZ(1,1,:);
if N2==1
	vz=ZZ(1,:);
end

[XX YY ZZ]=meshgrid(vx, vy,vz);
%XX=reshape(grid0,N0,N1,N2);
%YY=reshape(grid1,N0,N1,N2);
%ZZ=reshape(grid2,N0,N1,N2);

if which==0
return
end

if which==1 
	figure(10)
	hold;
	for i=1:size(pdiv, 1)-1
		st=pdiv(i)+1;
		en=pdiv(i+1);
		plot3(P0(st:en), P1(st:en), P2(st:en), 'b-');
	end
	plot3(grid0,grid1,grid2,'ko');
	xlabel('x axis (cm)');
	ylabel('y axis (cm)');
	zlabel('z axis (cm)');
	hold;
	return;
end


if which==2
	ar0=-r0/2:res0:r0/2-res0;
	ar1=-r1/2:res1:r1/2-res1;
	ar2=-r2/2:res2:r2/2-res2;
	
	vmax=round(max([max(B0(:)) max(B1(:)) max(B2(:))] ));
	vmin=round(min([min(B0(:)) min(B1(:)) min(B2(:))] ));
	
	
	v=vmin:10:vmax;

	figure(12)
	r2point=ceil(N2/2);
	r2pos=-r2/2+(r2point-1)*res2;
	[xx,yy]=meshgrid(ar0,ar1);
	subplot(3,3,1);
	[c,h]=contour(ar0,ar1,squeeze(B2(r2point,:,:)),v);
	%clabel(c,h);
	xlabel('x axis in cm');
	ylabel('y axis in cm');
	title(['Bz along the xy plane at z=' num2str(r2pos)]); 
	subplot(3,3,4); 
	[c,h]=contour(ar0,ar1,squeeze(B0(r2point,:,:)),v);
	%clabel(c,h);
	xlabel('r0 axis in cm');
	ylabel('r1 axis in cm');
	title(['Bx along the xy plane at z=' num2str(r2pos)]);
	subplot(3,3,7);
	[c,h]=contour(ar0,ar1,squeeze(B1(r2point,:,:)),v);
	%clabel(c,h);
	xlabel('r0 axis in cm');
	ylabel('r1 axis in cm');
	title(['By along the xy plane at z=' num2str(r2pos)]);


	r0point=ceil(N0/2);
	r0pos=-r0/2+(r0point-1)*res0;
	subplot(3,3,2)
	[c,h]=contour(ar1,ar2,squeeze(B2(:,:,r0point)),v);
	%clabel(c,h);
	xlabel('y axis in cm');
	ylabel('z axis in cm');
	title(['Bz along the yz plane at x=' num2str(r0pos)]);
	subplot(3,3,5)
	[c,h]=contour(ar1,ar2,squeeze(B0(:,:,r0point)),v);
	%clabel(c,h);
	xlabel('y axis in cm');
	ylabel('z axis in cm');
	title(['Bx along the yz plane at x=' num2str(r0pos)]);
	subplot(3,3,8)
	[c,h]=contour(ar1,ar2,squeeze(B1(:,:,r0point)),v);
	%clabel(c,h);
	xlabel('y axis in cm');
	ylabel('z axis in cm');
	title(['By along the yz plane at x=' num2str(r0pos)]);


	r1point=ceil(N1/2);
	r1pos=-r1/2+(r1point-1)*res1;
	subplot(3,3,3)
	[c,h]=contour(ar0,ar2,squeeze(B2(:,r1point,:)),v);
	%clabel(c,h);
	xlabel('x axis in cm');
	ylabel('z axis in cm');
	title(['Bz along the xz plane at y=' num2str(r1pos)]);
	subplot(3,3,6)
	[c,h]=contour(ar0,ar2,squeeze(B0(:,r1point,:)),v);
	%clabel(c,h);
	xlabel('x axis in cm');
	ylabel('z axis in cm');
	title(['Bx along the xz plane at y=' num2str(r1pos)]);
	subplot(3,3,9)
	[c,h]=contour(ar0,ar2,squeeze(B1(:,r1point,:)),v);
	%clabel(c,h);
	xlabel('x axis in cm');
	ylabel('z axis in cm');
	title(['By along the xz plane at y=' num2str(r1pos)]);


	br2p=uicontrol('Style','pushbutton','Units','normalized','Position',[.01 .66 .03 .03],'String','+');
	br2m=uicontrol('Style','pushbutton','Units','normalized','Position',[.05 .66 .03 .03],'String','-');
	br0p=uicontrol('Style','pushbutton','Units','normalized','Position',[.01 .36 .03 .03],'String','+');
	br0m=uicontrol('Style','pushbutton','Units','normalized','Position',[.05 .36 .03 .03],'String','-');
	br1p=uicontrol('Style','pushbutton','Units','normalized','Position',[.01 .07 .03 .03],'String','+');
	br1m=uicontrol('Style','pushbutton','Units','normalized','Position',[.05 .07 .03 .03],'String','-');

	%'quiver3(xx,yy,zeros(N0,N1),squeeze(B0(r2point,:,:)),squeeze(B1(r2point,:,:)),squeeze(B2(r2point,:,:)));' ...

	sr2=['r2pos=-r2/2+(r2point-1)*res2;' ...
	'subplot(3,3,1); ' ...
	'[c,h]=contour(ar0,ar1,squeeze(B2(r2point,:,:)),v);' ...
	'clabel(c,h);' ...
	'hold on;' ...
	'view(0,90);' ...
	'xlabel(''x axis in cm'');' ...
	'ylabel(''y axis in cm'');' ...
	'title([''Bz along the xy plane at z='' num2str(r2pos)]);' ...
	'hold off;' ...
	'subplot(3,3,4); ' ...
	'[c,h]=contour(ar0,ar1,squeeze(B0(r2point,:,:)),v);' ...
	'clabel(c,h);' ...
	'xlabel(''x axis in cm'');' ...
	'ylabel(''y axis in cm'');' ...
	'title([''Bz along the xy plane at z='' num2str(r2pos)]);' ...
	'subplot(3,3,7);' ...
	'[c,h]=contour(ar0,ar1,squeeze(B1(r2point,:,:)),v);' ...
	'clabel(c,h);' ...
	'xlabel(''x axis in cm'');' ...
	'ylabel(''y axis in cm'');' ...		
	'title([''Bz along the xy plane at z='' num2str(r2pos)]);'];

	sr0=['r0pos=-r0/2+(r0point-1)*res0;' ...
	'subplot(3,3,2);' ...
	'[c,h]=contour(ar1,ar2,squeeze(B2(:,:,r0point)),v);' ...
	'clabel(c,h);' ...
	'xlabel(''y axis in cm'');' ...
	'ylabel(''z axis in cm'');' ...
	'title([''Bz along the yz plane at x='' num2str(r0pos)]);' ...
	'subplot(3,3,5);' ...
	'[c,h]=contour(ar1,ar2,squeeze(B0(:,:,r0point)),v);' ...
	'clabel(c,h);' ...
	'xlabel(''y axis in cm'');' ...
	'ylabel(''z axis in cm'');' ...
	'title([''Bx along the yz plane at x='' num2str(r0pos)]);' ...
	'subplot(3,3,8);' ...
	'[c,h]=contour(ar1,ar2,squeeze(B1(:,:,r0point)),v);' ...
	'clabel(c,h);' ...
	'xlabel(''y axis in cm'');' ...
	'ylabel(''z axis in cm'');' ...
	'title([''Bx along the yz plane at x='' num2str(r0pos)]);'];

	sr1=['r1pos=-r1/2+(r1point-1)*res1;' ...
	'subplot(3,3,3);' ...
	'[c,h]=contour(ar0,ar2,squeeze(B2(:,r1point,:)),v);' ...
	'clabel(c,h);' ...
	'xlabel(''x axis in cm'');' ...
	'ylabel(''z axis in cm'');' ...
	'title([''Bz along the xz plane at r1='' num2str(r1pos)]);' ...
	'subplot(3,3,6);' ...
	'[c,h]=contour(ar0,ar2,squeeze(B0(:,r1point,:)),v);' ...
	'clabel(c,h);' ...
	'xlabel(''x axis in cm'');' ...
	'ylabel(''z axis in cm'');' ...
	'title([''Bx along the xz plane at r1='' num2str(r1pos)]);' ...
	'subplot(3,3,9);' ...
	'[c,h]=contour(ar0,ar2,squeeze(B1(:,r1point,:)),v);' ...
	'clabel(c,h);' ...
	'xlabel(''x axis in cm'');' ...
	'ylabel(''z axis in cm'');' ...
	'title([''By along the xz plane at r1='' num2str(r1pos)]);'];

	set(br2p,'Callback',['if (r2point<N2) r2point=r2point+1; end;',sr2]);
	set(br2m,'Callback',['if (r2point>1)  r2point=r2point-1; end;',sr2]);
	set(br0p,'Callback',['if (r0point<N0) r0point=r0point+1; end;',sr0]);
	set(br0m,'Callback',['if (r0point>1)  r0point=r0point-1; end;',sr0]);
	set(br1p,'Callback',['if (r1point<N1) r1point=r1point+1; end;',sr1]);
	set(br1m,'Callback',['if (r1point>1)  r1point=r1point-1; end;',sr1]);
	return;

elseif which==3
	
	
	mingz=min(min(min(ZZ)));
	maxgz=max(max(max(ZZ)));
	mingx=min(min(min(XX)));
	maxgx=max(max(max(XX)));
	mingy=min(min(min(YY)));
	maxgy=max(max(max(YY)));
	
	slices=5;
	lineslice=10;
	mview='a';
	
	if nargin>2
		slices=isoslices;
	end
	
	if nargin>3
		lineslice=numlines;
	end

	v=[-1000 -900 -700 -500 -300 -150 0 150 300 500 700 900 1000];
	
	if nargin>4
		mview=view;
	end


	figure(19);
	
	if mview=='a'
		subplot(1,3,1);
	end
	
	if mview=='a' | mview=='x'
		maxb0=max(max(max(B0)));
		minb0=min(min(min(B0)));

		v=minb0:abs(maxb0-minb0)/slices:maxb0;
		for i=1:size(v,2)-1
			p=patch(isosurface(XX,YY,ZZ,B0, v(i)));
			set(p, 'FaceColor', [0 i/size(v,2) 0], 'EdgeColor', 'none');
		end
		title(['Bx '])
		xlabel('x(cm)');
		ylabel('y(cm)'); 
		zlabel('z(cm)'); 

		[sx,sy,sz] = meshgrid(maxgx, mingy:(maxgy-mingy)/lineslice:maxgy,mingz:(maxgz-mingz)/lineslice:maxgz);
		hlines = streamline(XX,YY,ZZ,B0,B1,B2,sx,sy,sz);
		set(hlines,'LineWidth',2,'Color',[0.7 0 0]);

		[sx,sy,sz] = meshgrid(mingx, mingy:(maxgy-mingy)/lineslice:maxgy,mingz:(maxgz-mingz)/lineslice:maxgz);
		hlines = streamline(XX,YY,ZZ,B0,B1,B2,sx,sy,sz);
		set(hlines,'LineWidth',2,'Color',[0.7 0 0]);

		camlight;
		lighting gouraud;
		daspect([1, 1, 2]);
	end
	
	if mview=='a'
		subplot(1,3,2);
	end
	
	if mview=='a' | mview=='y'
		maxb1=max(max(max(B1)));
		minb1=min(min(min(B1)));
		v=minb1:abs(maxb1-minb1)/slices:maxb1;
		for i=1:size(v,2)-1
			p=patch(isosurface(XX,YY,ZZ,B1, v(i)));
			set(p, 'FaceColor', [0 0 i/size(v,2)], 'EdgeColor', 'none');
		end
		title(['By '])
		xlabel('x(cm)');
		ylabel('y(cm)'); 
		zlabel('z(cm)'); 

		[sx,sy,sz] = meshgrid(mingx:(maxgx-mingx)/lineslice:maxgx, maxgy,mingz:(maxgz-mingz)/lineslice:maxgz);
		hlines = streamline(XX,YY,ZZ,B0,B1,B2,sx,sy,sz);
		set(hlines,'LineWidth',2,'Color',[0.9 0.9 0]);

		[sx,sy,sz] = meshgrid(mingx:(maxgx-mingx)/lineslice:maxgx, mingy,mingz:(maxgz-mingz)/lineslice:maxgz);
		hlines = streamline(XX,YY,ZZ,B0,B1,B2,sx,sy,sz);
		set(hlines,'LineWidth',2,'Color',[0.9 0.9 0]);

		camlight;
		lighting gouraud;
		daspect([1, 1, 2]);
	end
	
	if mview=='a'
		subplot(1,3,3);
	end
	
	if mview=='a' | mview=='z'
		maxb2=max(max(max(B2)));
		minb2=min(min(min(B2)));


		v=minb2:abs(maxb2-minb2)/slices:maxb2;
		for i=1:size(v,2)-1
			p=patch(isosurface(XX,YY,ZZ,B2, v(i)));
			set(p, 'FaceColor', [i/size(v,2) 0 0], 'EdgeColor', 'none');

		end
		title(['Bz '])
		xlabel('x(cm)');
		ylabel('y(cm)'); 
		zlabel('z(cm)'); 

		[sx,sy,sz] = meshgrid(mingx:(maxgx-mingx)/lineslice:maxgx, mingy:(maxgy-mingy)/lineslice:maxgy,maxgz);
		hlines = streamline(XX,YY,ZZ,B0,B1,B2,sx,sy,sz);
		set(hlines,'LineWidth',2,'Color',[0 .7 0]);

		[sx,sy,sz] = meshgrid(mingx:(maxgx-mingx)/lineslice:maxgx, mingy:(maxgy-mingy)/lineslice:maxgy,mingz);
		hlines = streamline(XX,YY,ZZ,B0,B1,B2,sx,sy,sz);
		set(hlines,'LineWidth',2,'Color',[0 .7 0]);

		camlight;
		lighting gouraud;
		daspect([1, 1 ,2]);
	end
	return;
elseif which==4
	
	mingz=min(ZZ(:));
	maxgz=max(ZZ(:));
	mingx=min(XX(:));
	maxgx=max(XX(:));
	mingy=min(YY(:));
	maxgy=max(YY(:));

	slices=4;
	lineslice=10;
	
	plotc=1;
	mview='a';
	
	if nargin>2
		slices=isoslices;
	end
	
	if nargin>3
		lineslice=numlines;
	end

	if nargin>5
		plotc=plotcoil;
	end
	if nargin>4
		mview=view;
	end

	figure(100);
	
	v=(mingx+abs(maxgx-mingx)/slices):abs(maxgx-mingx)/(slices):(maxgx-abs(maxgx-mingx)/slices);
	if size(v,2)==0
		v=0;
	end
	
	if mview=='a'
		subplot(1,3,1);
	end
	
	if mview=='a' | mview=='x'
		hold on;	
		if size(findstr(vr,'R12'), 1) ~= 0
			h=slice(XX,YY,ZZ, B0, v, [],[]);
			set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8,'FaceAlpha', 0.9);
			hxm=slice(XX,YY,ZZ, B0, maxgx, [],[]);
			set(hxm,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8,'FaceAlpha', 0.5);
			hxmin=slice(XX,YY,ZZ, B0, mingx, [],[]);
			set(hxmin,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8,'FaceAlpha', 0.5);
			hyz=slice(XX,YY,ZZ, B0, [], maxgy,mingz);
			set(hyz,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.2);
		else
			h=slice(XX,YY,ZZ, B0, v, [],[]);
			set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
			hxm=slice(XX,YY,ZZ, B0, maxgx, [],[]);
			set(hxm,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
			hxmin=slice(XX,YY,ZZ, B0, mingx, [],[]);
			set(hxmin,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
			hyz=slice(XX,YY,ZZ, B0, [], maxgy,mingz);
			set(hyz,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.2);
		end


		hcont = ...
		contourslice(XX,YY,ZZ,B0,v,maxgy,mingz);
		set(hcont,'EdgeColor',[.4,.4,.4],'LineWidth',0.5);

		colormap jet

		[sx,sy,sz] = meshgrid(maxgx, mingy:(maxgy-mingy)/lineslice:maxgy,mingz:(maxgz-mingz)/lineslice:maxgz);
		hlines = streamline(XX,YY,ZZ,B0,B1,B2,sx,sy,sz);
		set(hlines,'LineWidth',2,'Color',[0.7 0 0]);

		[sx,sy,sz] = meshgrid(mingx, mingy:(maxgy-mingy)/lineslice:maxgy,mingz:(maxgz-mingz)/lineslice:maxgz);
		hlines = streamline(XX,YY,ZZ,B0,B1,B2,sx,sy,sz);
		set(hlines,'LineWidth',2,'Color',[0.7 0 0]);

		if plotc==1; 
			if size(pdiv, 1)~=1 & pdiv(1) ~= 0;
				for i=1:size(pdiv, 1)-1
					st=pdiv(i)+1;
					en=pdiv(i+1);
					p=line(P0(st:en), P1(st:en), P2(st:en), 'LineWidth', 3, 'Color', [0 0 0]);
				end
			end
		end
		colorbar;
		title(['Bx '])
		xlabel('x(cm)');
		ylabel('y(cm)'); 
		zlabel('z(cm)'); 
		%daspect([1,1,2])
		material shiny;
		axis tight
		box on;
		lighting phong;
		camlight(-37.5, 30, 'local');
		camlight(-37.5, 30, 'local');
	end
	
	if mview=='a'
		subplot(1,3,2);
	end
	
	if mview=='a' | mview=='y'
		hold on

		if size(findstr(vr,'R12'), 1) ~= 0
			h=slice(XX,YY,ZZ, B1, v, [],[]);
			set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8,'FaceAlpha', 0.9);
			hxm=slice(XX,YY,ZZ, B1, maxgx, [],[]);
			set(hxm,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8,'FaceAlpha', 0.5);
			hxmin=slice(XX,YY,ZZ, B1, mingx, [],[]);
			set(hxmin,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8,'FaceAlpha', 0.5);
			hyz=slice(XX,YY,ZZ, B1, [], maxgy,mingz);
			set(hyz,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.2);
		else
			h=slice(XX,YY,ZZ, B1, v, [],[]);
			set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
			hxm=slice(XX,YY,ZZ, B1, maxgx, [],[]);
			set(hxm,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
			hxmin=slice(XX,YY,ZZ, B1, mingx, [],[]);
			set(hxmin,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
			hyz=slice(XX,YY,ZZ, B1, [], maxgy,mingz);
			set(hyz,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.2);
		end

		hcont = contourslice(XX,YY,ZZ,B1,v,maxgy,mingz);
		set(hcont,'EdgeColor',[.4,.4,.4],'LineWidth',0.5);

		colormap jet

		[sx,sy,sz] = meshgrid(mingx:(maxgx-mingx)/lineslice:maxgx, maxgy,mingz:(maxgz-mingz)/lineslice:maxgz);
		hlines = streamline(XX,YY,ZZ,B0,B1,B2,sx,sy,sz);
		set(hlines,'LineWidth',2,'Color',[0.9 0 0]);

		[sx,sy,sz] = meshgrid(mingx:(maxgx-mingx)/lineslice:maxgx, mingy,mingz:(maxgz-mingz)/lineslice:maxgz);
		hlines = streamline(XX,YY,ZZ,B0,B1,B2,sx,sy,sz);
		set(hlines,'LineWidth',2,'Color',[0.9 0 0]);
		
		if plotc==1
			for i=1:size(pdiv, 1)-1
				st=pdiv(i)+1;
				en=pdiv(i+1);
				p=line(P0(st:en), P1(st:en), P2(st:en), 'LineWidth', 3, 'Color', [0 0 0]);
			end
		end

		colorbar;
		title(['By '])
		xlabel('x(cm)');
		ylabel('y(cm)'); 
		zlabel('z(cm)'); 
		%daspect([1,1,2])
		material shiny;
		axis tight
		box on;
		lighting phong;
		camlight(-37.5, 30, 'local');
		camlight(-37.5, 30, 'local');
	end
	
	if mview=='a'
		subplot(1,3,3);
	end
	
	if mview=='a' | mview=='z'
	
		hold on;
		%B2(1,1,1)=B2(1,1,1)-1e-10;
		if size(findstr(vr,'R12'), 1) ~= 0
			h=slice(XX,YY,ZZ, B2, v, [],[]);
			set(h,'FaceColor','interp','FaceLighting', 'flat','EdgeColor','none','DiffuseStrength',.8,'FaceAlpha', 0.9);
			hxm=slice(XX,YY,ZZ, B2, maxgx, [],[]);
			set(hxm,'FaceColor','interp','FaceLighting', 'flat','EdgeColor','none','DiffuseStrength',.8,'FaceAlpha', 0.5);
			hxmin=slice(XX,YY,ZZ, B2, mingx, [],[]);
			set(hxmin,'FaceColor','interp','FaceLighting', 'flat','EdgeColor','none','DiffuseStrength',.8,'FaceAlpha', 0.5);
			hyz=slice(XX,YY,ZZ, B2, [], maxgy,mingz);
			set(hyz,'FaceColor','interp','FaceLighting', 'flat','EdgeColor','none','DiffuseStrength',.2);
		else
			h=slice(XX,YY,ZZ, B2, v, [],[]);
			set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
			hxm=slice(XX,YY,ZZ, B2, maxgx, [],[]);
			set(hxm,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
			hxmin=slice(XX,YY,ZZ, B2, mingx, [],[]);
			set(hxmin,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
			hyz=slice(XX,YY,ZZ, B2, [], maxgy,mingz);
			set(hyz,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.2);
		end
		if lineslice>0
			[sx,sy,sz] = meshgrid(mingx:(maxgx-mingx)/lineslice:maxgx, mingy:(maxgy-mingy)/lineslice:maxgy,maxgz);
			hlines = streamline(XX,YY,ZZ,B0,B1,B2,sx,sy,sz);
			set(hlines,'LineWidth',2,'Color',[.9 0 0]);

			[sx,sy,sz] = meshgrid(mingx:(maxgx-mingx)/lineslice:maxgx, mingy:(maxgy-mingy)/lineslice:maxgy,mingz);
			hlines = streamline(XX,YY,ZZ,B0,B1,B2,sx,sy,sz);
			set(hlines,'LineWidth',2,'Color',[0.9 0 0]);
		end
	
		hcont = ...
		contourslice(XX,YY,ZZ,B2,v,maxgy,mingz);
		set(hcont,'EdgeColor',[.4,.4,.4],'LineWidth',0.5);

		colormap jet
		if plotc==1
			if size(pdiv, 1) ~=0 & pdiv(1)~=0
				for i=1:size(pdiv, 1)-1
					st=pdiv(i)+1;
					en=pdiv(i+1);
					p=line(P0(st:en), P1(st:en), P2(st:en), 'LineWidth', 3, 'Color', [0 0 0]);
				end
			end
		end
		colorbar;
		title(['Bz '])
		xlabel('x(cm)');
		ylabel('y(cm)'); 
		zlabel('z(cm)'); 
		%daspect([1,1,2])
		material shiny;
		axis tight
		box on;
		lighting phong;
		camlight(-37.5, 30, 'local');
		camlight(-37.5, 30, 'local');
	end
end

