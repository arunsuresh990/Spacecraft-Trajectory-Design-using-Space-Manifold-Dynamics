function dynamicinversionc2
% Problem parameters
ceta=.8;
ts=5;
w= 4/(ceta*ts);
%w= .9714;
y0(1)=0.1;
y0(2)=1;
y0(3)=0;
numSteps=1000;
t=linspace(0,10,numSteps)'
xref=(t);
plot(t(:,1),xref(:,1),'r');
options = odeset('RelTol',1e-12,'AbsTol',1e-8);             
[t,y] = ode45(@f,t,y0,options);
p=y
hold on
plot(t(:,1),y(:,1),'b');
grid on
title('Dynamic Inversion');
ylabel('x,xref');
xlabel('t');
figure
plot(t(:,1),y(:,3),'b');
function dydt = f(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
xr=t;
xr1=1;
xr2=0;
u=(1/y(1)+y(2)^2)*((xr2+y(1)*y(2)^2)-2*ceta*w*(y(2)-xr1)-w^2*(y(1)-xr));

dydt=[y(2);
     -y(1)*y(2)^2+(y(1)+y(2)^2)*u;
     u];
  end  


end  
  
  