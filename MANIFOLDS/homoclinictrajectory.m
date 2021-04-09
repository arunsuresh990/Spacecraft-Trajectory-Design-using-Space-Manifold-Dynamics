function homoclinictrajectory(i0)
format long;
mu =2.523e-5;
mu1=1-mu;
epsilon=2*10^(-5);

%i0=[  -0.850347046128398  -0.365889817560588  -0.010556064438523  -0.1635];
c= jacobi(i0(1),i0(2),i0(3),i0(4),mu)%3.00576
title(c);
xlabel(mu);
hillscurve(c);
tspan=[0 17.75];
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y1]=ode113(@f1,tspan,i0,options);
title(c);
xlabel(mu);
plot(y1(:,1),y1(:,2),'g');
hold on;
;


function dydt = f(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(2))^2);
r2=sqrt((-1+mu+y(1))^2+(y(2))^2);
dydt=[y(3);
      y(4);
      (y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      (y(2)-2*y(3)-(1-mu)*y(2)/(r1)^3-mu*y(2)/(r2)^3)];
  end  
function dydt = f1(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(2))^2);
r2=sqrt((-1+mu+y(1))^2+(y(2))^2);
dydt=[-y(3);
     - y(4);
      -(y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      -(y(2)-2*y(3)-(1-mu)*y(2)/(r1)^3-mu*y(2)/(r2)^3)];
end
end