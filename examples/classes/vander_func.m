

function dydt= WillRoss_func(t, iny)
	R=1000.0; 
	dydt=zeros(2,1);
	dydt(1)=iny(2);
	dydt(2)=R*(1.0-iny(1)*iny(1))*iny(2)-iny(1);
	
