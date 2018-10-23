Hello Bo,

Here are the beasts!!!

I define the spatial parts and the time dependences below!!

Co=1/8*sqrt(3)/sqrt(2)*(3*cos(thetap)*cos(thetap)-1+etap*sin(thetap)*sin(thetap)*cos(2*phip));
C1=1/4*sqrt(3/2)*sin(2*thetap)*(1-1/3*etap*cos(2*phip));
C2=1/8*sqrt(3/2)*(sin(thetap)*sin(thetap)+etap/3*(cos(thetap)*cos(thetap)+1)*cos(2*phip));
S1=1/2/sqrt(6)*etap*sin(thetap)*sin(2*phip);
S2=1/4/sqrt(6)*etap*cos(thetap)*sin(2*phip);

A20=-sqrt(3/2)*C1*2*sqrt(2)/3*Cos(wrt+phil)-sqrt(3/2)*2*sqrt(2)/3*S1*Sin(wrt+phil)+sqrt(6)C2*2/3*cos(2*wrt+2*phil)-sqrt(6)*2/3*sin(2*wrt+2*phil);

A21=Co*2*sqrt(2)/3+C1*(-1/3*cos(wrt+phil)-i/sqrt(3)*sin(wrt+phil))+S1*(-1/3*sin(wrt+phil)+i*1/sqrt(3)*cos(wrt+phil))-C2*(2*sqrt(2)/3*cos(2*wrt+2*phil)-2*i*sqrt(2/3)*sin(2*wrt+2*phil))+S2*(2*sqrt(2)/3*sin(2*wrt+2*phil)+2*i*sqrt(2/3)*cos(2*wrt+2*phil));

A21m=-Co*2*sqrt(2)/3-C1*(-1/3*cos(wrt+phil)+i/sqrt(3)*sin(wrt+phil))-S1*(-1/3*sin(wrt+phil)-i*1/sqrt(3)*cos(wrt+phil))+C2*(2*sqrt(2)/3*cos(2*wrt+2*phil)+2*i*sqrt(2/3)*sin(2*wrt+2*phil))-S2*(2*sqrt(2)/3*sin(2*wrt+2*phil)-2*i*sqrt(2/3)*cos(2*wrt+2*phil)

A22=Co*2/3+C1*(sqrt(2)/3*cos(wrt+phil)-i*sqrt(2/3)*sin(wrt+phil))+S1*(sqrt(2)/3*sin(wrt+phil)+i*sqrt(2/3)*cos(wrt+phil))+C2*(4/3*cos(2*wrt+2*phil)-2*i*1/sqrt(3)*sin(2*wrt+phil))-S2*(4/3*sin(2*wrt+phil)+2*i*1/sqrt(3)*cos(2*wrt+phil))


A22m=Co*2/3+C1*(sqrt(2)/3*cos(wrt+phil)+i*sqrt(2/3)*sin(wrt+phil))+S1*(sqrt(2)/3*sin(wrt+phil)-i*sqrt(2/3)*cos(wrt+phil))+C2*(4/3*cos(2*wrt+2*phil)+2*i*1/sqrt(3)*sin(2*wrt+phil))-S2*(4/3*sin(2*wrt+phil)-2*i*1/sqrt(3)*cos(2*wrt+phil))

H=wq*A20*(3*IZ*IZ)-3/2*wq*wq/wo*A21*A21m*[IZI- + I-IZ,IZI+ +I+IZ] - 3/2*wq*wq/wo*A22*A22m*[(I-)^2,(I+)^2]  -sqrt(3/2)*wq*wq/wo*3([IZ*IZ,IZI+ + I+IZ]A20*A21 - A21m*A20*[IZ*IZ,IZI- + I-IZ] -A20*A22/2*[IZ*IZ,(I+)^2] + A20*A22m/2*[IZ*IZ,(I-)^2])


Or in terms of Spherical tensors, If you call HQ= Sum_{m}A2mT2-m(-1)^m,

The hamiltonian I'd like to see is simply H=A20T20 - Sum_{m>0}[T2-m,T2m]/m*A2m*A2-m + Sum_{m ~= 0}[T20,T2m]/m*A20*A2m.


Thanks Bo!!!

