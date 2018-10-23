

ggg=0:9;
figure(10);
hold on;
ct=1;
loo=jet(length(ggg));
fbase='grid';
colordef black;
for gg=ggg
	cfN=fbase;
	cfN=strcat(cfN, num2str(gg));
	load(cfN);
	cfG=eval(cfN);
	h=plot3(cfG(:,1), cfG(:,2), cfG(:,3), 'o');
	set(h, 'Color', loo(ct,:));
	ct=ct+1;
end
hold off

	
