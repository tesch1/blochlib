
function jac=WillRoss_jacobian(t, iny)
	S=30.0; R=0.415; B=16.5; C=10.0;
	jac=zeros(3,3);
	jac(1,1)=S-R*iny(1)-iny(2)-iny(3);
	jac(1,2)=-iny(1);
	jac(1,3)=-iny(1);

	jac(2,1)=iny(2);
	jac(2,2)=iny(1)-C;
	jac(2,3)=0;

	jac(3,1)=-iny(3);
	jac(3,2)=0;
	jac(3,3)=B-iny(1)-iny(3);
	
