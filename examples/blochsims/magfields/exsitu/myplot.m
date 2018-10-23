
base='field1_';
numGrids=3;
start=0;
whic='x';
ct=1;
tit={'X-Y' 'X-Z' 'Y-Z'};
arr=['x' 'y' 'x' 'z' 'y' 'z'];
for i=start:(numGrids-1)
	[tx ty tz tbx tby tbz]=plotmag(strcat(base, num2str(i),'.mat'), 0);
	%if i==0
	%	tx=tx(1,:,:);
	%	ty=tz(1,:,:);
	%	tz=tz(1,:,:);
	%	tbx=tbx(1,:,:);
	%	tby=tby(1,:,:);
	%	tbz=tbz(1,:,:);
	%end
	eval(strcat('x', num2str(i), '=squeeze(tx);'));
	eval(strcat('y', num2str(i), '=squeeze(ty);'));
	eval(strcat('z', num2str(i), '=squeeze(tz);'));
	%tmax=max(max(max(abs(tbx(:))), max(abs(tby(:)))), max(abs(tbz(:))));
	tmax=1;
	eval(strcat('bx', num2str(i), '=squeeze(tbx)/tmax;'));
	eval(strcat('by', num2str(i), '=squeeze(tby)/tmax;'));
	eval(strcat('bz', num2str(i), '=squeeze(tbz)/tmax;'));
	
end

figure(15);
for i=start:(numGrids-1)
	crB=squeeze(eval(strcat('b',whic,num2str(i))));
	crX=eval(strcat(arr(i*2+1),num2str(i)));
	crY=eval(strcat(arr(i*2+2),num2str(i)));
	[tx ty]=gradient(crB,abs(crX(1,1)-crX(2,1)),abs(crY(1,1)-crY(1,2))); 
	eval(strcat('gx', num2str(i), '=tx;'));
	eval(strcat('gy', num2str(i), '=ty;'));
	
	crGx=eval(strcat('gx', num2str(i)));
	crGy=eval(strcat('gy', num2str(i)));
	
	subplot(3,3,ct);
	pcolor(crX,crY,crB); shading interp;  colorbar;
	title(strcat(tit(i+1), 'plane B', whic));
	xlabel(strcat(arr(i*2+1), ' (cm)'));
	ylabel(strcat(arr(i*2+2), ' (cm)'));
	
	subplot(3,3,ct+1);
	pcolor(crX,crY,crGx); shading interp;   colorbar;
	title(strcat(tit(i+1), 'plane d(B', whic,')/dx'));
	xlabel(strcat(arr(i*2+1), ' (cm)'));
	ylabel(strcat(arr(i*2+2), ' (cm)'));

	subplot(3,3,ct+2);
	pcolor(crX,crY,crGy); shading interp;   colorbar;
	title(strcat(tit(i+1), 'plane d(B', whic,')/dy'));
	xlabel(strcat(arr(i*2+1), ' (cm)'));
	ylabel(strcat(arr(i*2+2), ' (cm)'));

	ct=ct+3;
end
	

figure(21);
subplot(1,3,1);
hold on;
crB=squeeze(eval(strcat('b',whic,'0')));	
h=surf(x0,y0,z0,crB); shading interp;
set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
crB=squeeze(eval(strcat('b',whic,'1')));	
h=surf(x1,y1,z1,crB);shading interp;
set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
crB=squeeze(eval(strcat('b',whic,'2')));	
h=surf(x2,y2,z2,crB);shading interp;
set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
title(strcat('B',whic,' on slices'));
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
colorbar
hold off;

subplot(1,3,2);
hold on;
h=surf(x0,y0,z0,gx0); shading interp;
set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
h=surf(x1,y1,z1,gx1);shading interp;
set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
h=surf(x2,y2,z2,gx2);shading interp;
set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
title(strcat('d(B',whic,')/dx on slices'));
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
colorbar
hold off;

subplot(1,3,3);
hold on;
h=surf(x0,y0,z0,gy0); shading interp;
set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
h=surf(x1,y1,z1,gy1);shading interp;
set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
h=surf(x2,y2,z2,gy2);shading interp;
set(h,'FaceColor','interp','EdgeColor','none','DiffuseStrength',.8);
title(strcat('d(B',whic,')/dy on slices'));
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
colorbar
hold off;
