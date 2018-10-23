

function dydt= WillRoss_func(t, iny)
	S=30.0; R=0.415; B=16.5; C=10.0;
%dydt[0].x()=S*iny[0].x()-R*iny[0].x()*iny[0].x()-iny[0].x()*iny[0].y()-iny[0].x()*iny[0].z();
%dydt[0].y()=iny[0].x()*iny[0].y()-C*iny[0].y();
%dydt[0].z()=B*iny[0].z()-iny[0].x()*iny[0].z()-0.5*iny[0].z()*iny[0].z();

	dydt=zeros(3,1);
	dydt(1)=S*iny(1)-R*iny(1)*iny(1)-iny(1)*iny(2)-iny(1)*iny(3);
	dydt(2)=iny(1)*iny(2)-C*iny(2);
	dydt(3)=B*iny(3)-iny(1)*iny(3)-0.5*iny(3)*iny(3);
	
