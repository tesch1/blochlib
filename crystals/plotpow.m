load rep256;
dd(:,3)=1;

[X,Y,Z]=sph2cart(dd(:,1), dd(:,2),1);

figure(1);
plot3(X,Y,Z, 'ko');


figure(2);
plot3(cos(dd(:,1)).*cos(dd(:,2)), ...
	sin(dd(:,1)).*cos(dd(:,2)), ...
	sin(dd(:,2)), 'ko');




