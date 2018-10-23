function h=plottrag(data,type, shows, vmat)

% function h=plottrag(data,type, shows, vmat)
%
% author:: Bo Blanton UC Berkeley, Dept of Chemistry
% email:: magneto@dirac.cchem.berkeley.edu
% more info:: http://waugh.cchem.berkeley.edu/blochlib/
% last modified:: 10.20.02
%
% plots and returns the data from a trajectory file
% from the 'blochlib' class "BlochSolver.writeData(ofstream)"
% or from a continous write process "BlochSolver.setWritePolicy(SolverOps::Continous)"
%
% the data should be organized like
% <spin1 X> <spin1 Y> <spin1 Z> ... <spinN X> <spinN Y> <spinN Z>
%
% data-->the input data
% type--> type of plot...either 'BlochSphere[1]' or 'time v Mx v My'[2]
% shows--> if a single number (i.e. 3) shows that spins trajectory
%	   if a vector (i.e. 1:30) shows those spins...
%           if 'a' shows all the trajectory
%           if not present will show all the trajectories
% 
% vmat--> AN ARRAY...iff not preset will display every point
%         inputting "1:50" will show the first 50 points of tragectory

if nargin<2; type=1; end;
if nargin<3; shows='a'; end;
if nargin<4; vmat=1:size(data, 1); end;

if vmat(1) < 0
	disp('input view size has start less then 1, fixing');
	vmat=1:vmat(length(vmat));
end

if length(vmat) > size(data, 1)
	disp('input view size is too long for input data, fixing');
	vmat=vmat(1):size(data, 1);
end

if vmat(length(vmat)) > size(data, 1)
	disp('input view size is too long for input data, fixing');
	vmat=vmat(1):size(data, 1);
end

h=figure(30);

if type==1
	clf reset;
	colordef white;
	maxV=max(max(max(abs(data))));
	spsize=39;
	stp=10;
	[xc,yc, zc]=sphere(spsize);
	
	
	xc=xc*maxV;
	yc=yc*maxV;
	zc=zc*maxV;
	
	fcolor=[0.8, 0, 0];
	lcolor1=[0.50, 0.50, 0.7];
	lcolor2=[0.50, 0.50, 0.7];
	%xcc=[xc(:,1); xc(:,20)]; ycc= [yc(:,1); yc(:,20)]; zcc=[zc(:,1); zc(:,20)];
	%line(xcc,ycc,zcc ,'Color',fcolor);
	ff=fill3(xc(20,:), yc(20,:), zc(20,:),bone(size(xc(20,:),1)));
	%set(ff, 'Cdata', jet(size(xc(20,:))));
	%set(ff, 'FaceAlpha', 0.4);
	ff=line(xc(:,1:stp:((spsize+1)/2)),yc(:,1:stp:((spsize+1)/2)),zc(:,1:stp:((spsize+1)/2)), 'Color', lcolor1);
	set(ff, 'LineStyle', ':');   
	ff=line(xc(:,((spsize+1)/2):stp:spsize),yc(:,((spsize+1)/2):stp:spsize),zc(:,((spsize+1)/2):stp:spsize), 'Color', lcolor1);
	xc=xc'; yc=yc'; zc=zc';
	%line(xc(:,1:stp:spsize),yc(:,1:stp:spsize),zc(:,1:stp:spsize), 'Color', lcolor2);
	
	if shows=='a';
		cols=bone(ceil(size(data, 2))/2);
		for i=0:ceil(size(data, 2)/3)-1
			line(data(vmat,(i*3)+1), data(vmat, (i*3)+2), data(vmat, (i*3)+3), 'Color', cols(i+1,:));
		end
	elseif max(size(shows))>2;
		for i=1:length(shows)
			line(data(vmat,((shows(i)-1)*3)+1), data(vmat, ((shows(i)-1)*3)+2), data(vmat, ((shows(i)-1)*3)+3));
		end
	else
		line(data(vmat,((shows-1)*3)+1), data(vmat, ((shows-1)*3)+2), data(vmat, ((shows-1)*3)+3));
	end	
	axis tight;
	view([1 1 1]);
	daspect([1 1 1]);
	box on;
	%set(h, 'Projection', 'perspective');
elseif type==2
	t=vmat;
	clf reset;
	colordef black;
	h=newplot;
	if shows=='a'
		loo=hot(round(size(data, 2)/3));
		for i=0:round(size(data, 2)/3)-1
			line(t, data(vmat,(i*3)+1), data(vmat, (i*3)+2), 'Color', loo(i+1,:));
		end
		view([1 1 1]);
	elseif length(shows)>1
		loo=jet(length(shows));
		for i=1:length(shows)
			line(t, data(vmat, ((shows(i)-1)*3)+1), data(vmat, ((shows(i)-1)*3)+2), 'Color', loo(i+1,:));
		end
		view([1 -1 1]);
	else
		line(t, data(vmat,((shows-1)*3)+1), data(vmat, ((shows-1)*3)+2));
	end
	axis tight;
	view([14 29]);
	%daspect([1 1 1]);
	box on;
end

