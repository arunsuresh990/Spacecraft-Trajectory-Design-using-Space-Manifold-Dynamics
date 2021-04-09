function dynamicinversion
% Problem parameters
ceta=.8;
ts=10;
w= 4/(ceta*ts);
y0(1)=1;
y0(2)=0;

tspan = [0 5];

options = odeset('RelTol',1e-12,'AbsTol',1e-8);             
[t,y] = ode45(@f,tspan,y0,options);
p=y
plot(t(:,1),y(:,1),'r');
hold on;
plot(t(:,1),y(:,2),'b');
grid on
title('Dynamic Inversion');
ylabel('x');
xlabel('t');

function dydt = f(t,y)
  % Derivative function -- mu and mustar shared with the outer function.

u=(-y(1)+y(1)^2)/(1+y(1));
dydt=[-y(1)^2+(1+y(1))*u;
     u];
  end  

  % -----------------------------------------------------------------------

 

  % -----------------------------------------------------------------------

end  % liaorbitode

  
  