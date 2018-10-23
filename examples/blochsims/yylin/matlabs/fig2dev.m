function Mp = fig2dev(t,M,gammah,td,tr,z)
nz=length(z);
M=reshape(M,nz,3);
Mp=zeros(size(M));
% -- DDF
B=[-(M(:,1)-mean(M(:,1))),-(M(:,2)-mean(M(:,2))),2*(M(:,3)-mean(M(:,3)))]/3/td/gammah;
% -- RD
B(:,1)=B(:,1)-mean(M(:,2))/tr(1)/gammah(1);
B(:,2)=B(:,2)+mean(M(:,1))/tr(1)/gammah(1);
% -- precession
Mp=gammah*cross(M,B);
% --- output
Mp=Mp(:);
M=M(:);
return
