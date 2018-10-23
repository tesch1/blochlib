function [fidt,fid,Y]=fig2gen(t,M,gammah,td,tr,z,prec)
i=sqrt(-1);
options=odeset('RelTol',prec,...
            'Stats','on');
M=M(:);
[T,Y]=ode45('fig2ode',t,M,options,gammah,td,tr,z);
fidt=T;
tmp=mean(reshape(Y,length(T),length(z),3),2);
fid=reshape(tmp(:),length(T),3);
return
