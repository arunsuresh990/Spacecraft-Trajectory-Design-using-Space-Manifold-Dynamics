function simpletrajectory
% Problem parameters
mu=1.215e-2;
mu1=1-mu;
c= 3.150;
hillscurve(c);
% y0(1)=1.283991897049027;
% y0(3)=-0.539090042308706;
% y0(2)=-1.354570087791366;
y0=[1.118284730835350 ; 0 ; 0; 0.185998679020517];
r1=sqrt((mu+y0(1))^2+(y0(2))^2);
r2=sqrt((-mu1+y0(1))^2+(y0(2))^2);
% y0(4)=-sqrt((y0(1)^2)+(y0(2)^2)+(2*(mu1)/r1)+(2*mu/r2)-(y0(3)^2+c));
% y0 = [y0(1);y0(2);y0(3); y0(4)]
hp=3.420555832978371;
tspan = [0 hp];

options = odeset('RelTol',1e-9,'AbsTol',1e-8);             
[t,y] = ode45(@f,tspan,y0,options);
p=y;
plot(y(:,1),y(:,2),'r');
grid on;
t1=t+hp;
t=[t;t+hp;t+2*hp;t+3*hp;t+4*hp]
y=[y;y;y;y;y];
figure
hillscurve(c);
plot(y(:,1),y(:,2),'b');
%axis([.9 1.05 -.03 .03]); 
title('Restricted three body problem');
ylabel('y(t)');
xlabel('x(t)');
hold on;
%plot (1.029,0,'*');
M=(1-mu);
J=-mu;
plot(J,0,'k*');
hold on;
plot(M,0,'b*');

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
 function dydt = f1(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(2))^2);
r2=sqrt((-1+mu+y(1))^2+(y(2))^2);
dydt=[-y(3);
     - y(4);
      -(y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      -(y(2)-2*y(3)-(1-mu)*y(2)/(r1)^3-mu*y(2)/(r2)^3)];
end
  % -----------------------------------------------------------------------

 

  % -----------------------------------------------------------------------

end  % liaorbitode

  
  