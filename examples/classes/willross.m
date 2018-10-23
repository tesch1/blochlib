

% a tester for the ode15s


tspan=0:0.01:50;

y0=[10 5 6];

opts=odeset('Jacobian', @WillRoss_jacobian, ...
            'Stats','on',...
            'JConstant','off',...
            'BDF', 'off');
tt=cputime;
figure(3);
[t, Y]=ode15s(@WillRoss_func, tspan, y0, opts);
plot3(Y(:,1), Y(:,2), Y(:,3),'b-');
box;
axis tight;
disp(cputime-tt)
