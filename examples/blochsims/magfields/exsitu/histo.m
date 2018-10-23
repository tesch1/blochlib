%the 'histogram plot  
load fieldHist.mat;
%subplot(3,1,1);
figure(33);
plot((B02-mean(B02))/mean(B02), B10, 'ko');
ylabel('B0z');
xlabel('B1x');
