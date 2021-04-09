function hetroclinictrajectory
format long;
mu=1.215e-2;
mu1=1-mu;
epsilon=2*10^(-4);

haloPeriod=2.844817312513691;
i0=[0.815962663523757;0;0; 0.207251594775628];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
k=40;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);

options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t1,y1]=ode113(@f,tspan,i0,options);%FUNCTION CALL
plot(y1(:,1),y1(:,2),'k');
hold on;

haloPeriod=3.420603969704802;
i0=[1.118284744395231;0;0;0.185998652126393];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
k=40;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,i0,options);%FUNCTION CALL
plot(y(:,1),y(:,2),'k');

hold on

% i0=[1-mu;-.07831;0.334381593505635;-.012];
% %i0=[1-mu;-.04971;0.538654291120618;-.0164];
% c= jacobi(i0(1),i0(2),i0(3),i0(4),mu)
% hillscurve(c);
% tspan=[0 2.8];
% options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
% [t,y1]=ode113(@f,tspan,i0,options);
% plot(y1(:,1),y1(:,2),'b');
% hold on;
% tspan=[0 1.69];
% options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
% [t,y]=ode113(@f1,tspan,i0,options)
% plot(y(:,1),y(:,2),'b');
% hold on;
% grid on;


i0=[0.814562011311580  ; 0.000119000140011 ;  0.014700952438085 ;  0.208213088752462];
c= jacobi(i0(1),i0(2),i0(3),i0(4),mu)
tspan=[0 3.4];
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y1]=ode113(@f,tspan,i0,options);
plot(y1(:,1),y1(:,2),'b');

hold on

i0=[0.813031511223778 ;  0.000257197053506  ; 0.013892141509814  ; 0.209974929051467];
c= jacobi(i0(1),i0(2),i0(3),i0(4),mu)
tspan=[0 3.8];
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y1]=ode113(@f,tspan,i0,options);
plot(y1(:,1),y1(:,2),'r');
grid on
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