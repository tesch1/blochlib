

minv=0.005;
r=0.03:0.0005:max(max(real(fid)));

ysw=(1:size(fid,1))*(0.000125);
swe=64; 
xsw=-swe/2:(swe/(size(fid,2)-1)):swe/2; 
[x,y]=meshgrid(xsw, ysw);

foo=zeros(size(fid));
for i=1:size(fid,1)
	for j=1:size(fid,2)
		if fid(i,j)>minv
			foo(i,j)=real(fid(i,j));
		end
	end
end

%pcolor(x,y,foo); shading interp;

contour(x,y,foo,r); axis([4 6 -Inf Inf])
