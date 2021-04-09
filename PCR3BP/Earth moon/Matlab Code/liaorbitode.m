function liaorbitode
% Finds an initial trajectory which could be used to generate a periodic
% orbit and hence all manifolds.This initial condition is passed to another
% program (liaorbitwwithnr(c,init1);) to obtain periodic orbit. Once
% periodic orbit is generated initial condition and time period is passed
% to next program to generate manifolds.
format long;
mu=1.215e-2;
%mu=0.012155
mu1=1-mu;
c= 3.1850;%3.00576%3.0075940671432188;
y0(1)= 0.829874;%00.81597%0.969539872159348, 1.0201059805843
y0(3)=0;
y0(2)=0;
r1=sqrt((mu+y0(1))^2+(y0(2))^2);
r2=sqrt((-mu1+y0(1))^2+(y0(2))^2);
y0(4)=sqrt((y0(1)^2)+(y0(2)^2)+(2*(mu1)/r1)+(2*mu/r2)-(y0(3)^2+c));%0.007530309627883
tspan = [0 40];
options = odeset('RelTol',1e-12,'AbsTol',1e-8,'Events',@events);             
[t,y,te,ye,ie] = ode45(@f,tspan,y0,options);
%plot(y(:,1),y(:,2),'r');
init1=y0(1);
%figure
liaorbitwwithnr(c,init1);
hold on;

plot (0.836918007316934,0,'*');
%plot (1.155679913094732,0,'*');
  
%hillscurve(c);
function dydt = f(t,y)
 
r1=sqrt((mu+y(1))^2+(y(2))^2);
r2=sqrt((-1+mu+y(1))^2+(y(2))^2);
dydt=[y(3);
      y(4);
      (y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      (y(2)-2*y(3)-(1-mu)*y(2)/(r1)^3-mu*y(2)/(r2)^3)];
  end  

  
  function [value,isterminal,direction] = events(t,y)
   value = y(2);
    isterminal = 1;        % stop at local minimum
    direction  = 1;         
  end


end  % liaorbitode

  
  