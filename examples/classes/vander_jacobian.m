
function jac=WillRoss_jacobian(t, iny)
	R=1000;
	jac=zeros(2,2);
	jac(1,1)=0.0;
	jac(1,2)=1.0;

	jac(2,1)=R*2.0*iny(1)*iny(2)-1.0;
	jac(2,2)=R*(1.0-iny(1)*iny(1));
	
