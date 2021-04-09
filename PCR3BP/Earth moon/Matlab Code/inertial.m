function inertial
% Problem parameters
c=input('Enter value of Energy,limits are "3.007640   3.007536   3.0000780   2.999921   2.999921":');
%hillscurve(c);
%grid on;
mu=1.215e-2;
mu1=1-mu;
plot (mu,0,'*');
%plot (mu1,0,'*');
hold on;
y0(1)=.285;
y0(2)=-.5050514;
y0(3)=-.010;
r1=sqrt((mu+y0(1))^2+(y0(2))^2);
r2=sqrt((-mu1+y0(1))^2+(y0(2))^2);
y0(4)=sqrt((y0(1)^2)+(y0(2)^2)+(2*(mu1)/r1)+(2*mu/r2)-c);%0.007530309627883
%ydot=sqrt(x^2+y^2+2*(1-mu)/r1 + 2*mu/r2-(xdot^2+C))
y0 = [y0(1);y0(2);y0(3); y0(4)];
tspan = [0 5];

options = odeset('RelTol',1e-12,'AbsTol',1e-8);             
[t,y] = ode45(@f,tspan,y0,options);
%plot(y(:,1),y(:,2),'r');
[m,n]=size(t);
p=t;
l=p*27.3*24*60*60/(2*pi);
omega=1.02e-5;
hold on;
for i=1:m
    theeta=omega*l(i,1)
    ix(i)=y(i,1)*cos(theeta)-y(i,2)*sin(theeta);
    iy(i)=y(i,1)*sin(theeta)+y(i,2)*cos(theeta);


plot(ix,iy,'b');
title('Restricted three body problem');
ylabel('y(t)');
xlabel('x(t)');
end

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
