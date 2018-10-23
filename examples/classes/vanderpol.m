

% a tester for the ode15s


tspan=0:2.5:10000;

y0=[2 0];

opts=odeset('Jacobian', @vander_jacobian, ...
            'OutputFcn',@odeplot,...
            'Stats','on',...
            'JConstant','off',...
            'BDF', 'on');
ton=cputime;
figure(3);
[t, Y]=ode15s(@vander_func, tspan, y0, opts);
disp(cputime-ton);
