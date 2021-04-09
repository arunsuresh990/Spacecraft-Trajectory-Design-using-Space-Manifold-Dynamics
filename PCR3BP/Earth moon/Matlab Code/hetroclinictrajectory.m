function hetroclinictrajectory
format long;
mu=1.215e-2;
mu1=1-mu;
epsilon=2*10^(-4);

haloPeriod= 3.079930437000846;
i0=[ 0.988717473332647;0;0; 0.010244790119764];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
k=40;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);

options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t1,y1]=ode113(@f,tspan,i0,options);%FUNCTION CALL
plot(y1(:,1),y1(:,2),'k');
hold on;

haloPeriod=3.117288267252508;
i0=[1.008126782443260;0;0;0.011128601253921];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
k=40;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,i0,options);%FUNCTION CALL
plot(y(:,1),y(:,2),'k');

hold on

i0=[mu1;6.33e-3;-0.011881326221789;0];
c= jacobimu(i0(1),i0(2),i0(3),i0(4),mu)
%hillscurve(c);
tspan=[0 3.2];
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y1]=ode113(@f,tspan,i0,options);
plot(y1(:,1),y1(:,2),'r');
hold on;
%tspan=[0 4.5];
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f1,tspan,i0,options);
plot(y(:,1),y(:,2),'r');
hold on;
plot (0.990027676383308,0,'*');
hold on;
plot (1.010033021348862,0,'*');
hold on;
plot (mu1,0,'*');




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