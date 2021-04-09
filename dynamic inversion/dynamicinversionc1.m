function dynamicinversionc1
% Problem parameters
ceta=.3;
ts=10;
w= 4/(ceta*ts);
%w= .9714;
y0(1)=500;
y0(2)=0;
y0(3)=0;

numSteps=1000;
i=numSteps;
t=linspace(0,25,numSteps)';
xref=(5+2*t+3*t.^2);
plot(t(:,1),xref(:,1),'k');
options = odeset('RelTol',1e-12,'AbsTol',1e-8);             
[t,y] = ode45(@f,t,y0,options);
y01=[10 2];
t1=linspace(0,25,numSteps)';
[t1,y1] = ode45(@f1,t,y01,options);
p=y
hold on
plot(t(:,1),y(:,1),'k:');
grid on
title('Dynamic Inversion');
ylabel('x,xref');
xlabel('t');
% figure
% plot(t(:,1),y(:,3),'b');
% 
% figure
% plot(t1(:,1),y1(:,1),'k');

function dydt = f(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
xr=5+2*t+3*t^2;
xr1=2+6*t;
xr2=6;
u=(1/2)*((xr2+y(1)^3)-2*ceta*w*(y(2)-xr1)-w^2*(y(1)-xr));
dydt=[y(2);
     -y(1)^3+2*u;
     u];

  end  
function dydt = f1(t,y)
  % Derivative function -- mu and mustar shared with the outer function.

dydt=[y(2);
     -y(1)^3+2];

  end  

end  
  
  