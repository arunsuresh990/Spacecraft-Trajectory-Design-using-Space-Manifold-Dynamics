function insidehills
% Problem parameters
c=input('Enter value of Energy,limits are "  3.188335717526626   3.172155838875999   3.012146565419430   2.987997622500000   2.987997622500000 ":');
hillscurve(c);
grid on;
mu=1.215e-2;
mu1=1-mu;
y0(1)=.7958;
y0(3)=.05050514;
y0(2)=-.0010;
r1=sqrt((mu+y0(1))^2+(y0(2))^2);
r2=sqrt((-mu1+y0(1))^2+(y0(2))^2);
y0(4)=sqrt((y0(1)^2)+(y0(2)^2)+(2*(mu1)/r1)+(2*mu/r2)-(y0(3)^2+c));%0.007530309627883
y0 = [y0(1);y0(2);y0(3); y0(4)]
tspan = [0 10];

options = odeset('RelTol',1e-12,'AbsTol',1e-8);             
[t,y] = ode45(@f,tspan,y0,options);
plot(y(:,1),y(:,2),'b');
title('Restricted three body problem');
ylabel('y(t)');
xlabel('x(t)');

  % -----------------------------------------------------------------------
  % Nested functions -- problem parameters provided by the outer function.
  %

  function dydt = f(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(2))^2);
r2=sqrt((-1+mu+y(1))^2+(y(2))^2);
dydt=[y(3);
      y(4);
      (y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      (y(2)-2*y(3)-(1-mu)*y(2)/(r1)^3-mu*y(2)/(r2)^3)];
  end  

  % -----------------------------------------------------------------------

 

  % -----------------------------------------------------------------------

end  % liaorbitode
