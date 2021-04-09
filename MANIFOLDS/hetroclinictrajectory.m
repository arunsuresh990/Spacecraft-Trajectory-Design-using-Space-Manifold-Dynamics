function hetroclinictrajectory
format long;
mu = 7.802e-5;
mu1=1-mu;
epsilon=2*10^(-4);

haloPeriod= 3.025210183985142;
i0=[ 0.967140945944572;0;0; 0.027284635429929];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
k=40;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);

options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t1,y1]=ode113(@f,tspan,i0,options);%FUNCTION CALL
plot(y1(:,1),y1(:,2),'k');
hold on;

haloPeriod=3.137767830378976;
i0=[ 1.025090755289090;0;0;0.027511316423054];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
k=40;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,i0,options);%FUNCTION CALL
plot(y(:,1),y(:,2),'k');

hold on

i0=[ 0.99992198;-.0195;0.018141538964882;-.019];
c= jacobi(i0(1),i0(2),i0(3),i0(4),mu)
hillscurve(c);
tspan=[0 1.5];
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y1]=ode113(@f,tspan,i0,options);
plot(y1(:,1),y1(:,2),'r');
hold on;
tspan=[0 1.2];
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f1,tspan,i0,options);
plot(y(:,1),y(:,2),'r');
hold on;


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