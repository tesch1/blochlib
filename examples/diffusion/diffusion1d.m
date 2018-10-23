
function diffusion1d(lenin, maxtin, iterin)

len=lenin;
maxx=1;
dx=maxx/len;

x=0:dx:maxx-dx;
x=x';
dxsq=dx*dx;

iter=iterin;
maxt=maxtin;
dt=maxt/iter;

diff=1;

fact=dt*diff/dxsq;

A=zeros(len,len);
B=zeros(len,len);

for i=1:len
	if i==1
		%A(i,i)=1;
		A(i,i)=1+fact;
		A(i, i+1)=-fact/2;
		
		%B(i,i)=1;
		B(i,i)=1-fact;
		B(i,i+1)=fact/2;
	elseif i==len
		
		%A(i,i)=1;
		A(i,i)=1+fact;
		A(i, i-1)=-fact/2;
		
		%B(i,i)=1;
		B(i,i)=1-fact;
		B(i,i-1)=fact/2;
	else
		A(i,i)=1+fact;
		A(i, i+1)=-fact/2;
		A(i, i-1)=-fact/2;
		
		B(i,i)=1-fact;
		B(i,i+1)=fact/2;
		B(i,i-1)=fact/2;
	end
end

data=zeros(iter+1, len);

ct=2;

sol=x;
data(1,:)=sol';
BCm1=x(1);
BC1=x(len);
B
dirch=0;
[L,U]=lu(A);
	
while ct<=iter+1
	sol=B*x;
	sol=U\(L\sol);
	x=sol;
	if dirch ==0
		 x(1)=BCm1; x(len)=BC1;	
	else
		x(1)=x(2); x(len)=x(len-1);	
	end
	data(ct,:)=x';
	ct=ct+1;
end

figure(1); pcolor(data'); shading interp;
figure(2);hold; plot(data(1,:), 'k');  plot(data(iter, :)); hold;


