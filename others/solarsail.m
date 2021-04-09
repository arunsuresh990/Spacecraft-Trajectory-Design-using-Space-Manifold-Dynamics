function solarsail
% Problem parameters
mu = 132712440018;

y0(1)= 1.49600000e8;
y0(2)= 0;
y0(3)=0 ;
%y0(4)=-0.043013862776183 ;

r=sqrt((y0(1))^2+(y0(2))^2);
y0(4)=sqrt(mu/r);
y0 = [y0(1);y0(2);y0(3); y0(4)];

tspan = [0 365*24*60*60*3];

options = odeset('RelTol',1e-12,'AbsTol',1e-8);             
[t,y] = ode45(@f,tspan,y0,options)
p=y;
plot(0,0,'b*');
hold on;
plot(1.49600000e8,0,'r*');
hold on;
plot(y(:,1),y(:,2),'b');


  function dydt = f(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
r=sqrt((y(1))^2+(y(2))^2);
dydt=[y(3);
      y(4);
      (-(mu/r^3)*y(1)-0.00000001356155);
      (-(mu/r^3)*y(2)-0.00000001356155)];
  end  


end  % liaorbitode

  
  