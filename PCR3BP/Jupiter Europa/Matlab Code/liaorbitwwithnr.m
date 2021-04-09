function liaorbitwwithnr(C,init)
%This program is mainly called bo Liaorbitode.
%To generate periodic orbit and manifolds (statements 50 to 56).
mu =2.523e-5;
mu1=1-mu;
delx=1e-10;
x1(1)=init;
for i=1:10
c=C;
y0(1)=x1(i);
y0(2)=0;
y0(3)=0;
r1=sqrt((mu+y0(1))^2+(y0(2))^2);
r2=sqrt((-mu1+y0(1))^2+(y0(2))^2);
y0(4)=sqrt((y0(1)^2)+(y0(2)^2)+(2*(mu1)/r1)+(2*mu/r2)-(y0(3)^2+c));%0.007530309627883
y0 = [y0(1);y0(2);y0(3); y0(4)];
tspan = [0 20];
options = odeset('RelTol',1e-11,'AbsTol',1e-22,'Events',@events);
             
[t,y,te,ye,ie] = ode45(@f,tspan,y0,options);
[m,n]=size(y);
xf1=y(m,1);
f1=x1(i)-xf1;
y0(1)=x1(i)+delx;
y0(2)=0;
y0(3)=0;
r1=sqrt((mu+y0(1))^2+(y0(2))^2);
r2=sqrt((-mu1+y0(1))^2+(y0(2))^2);
y0(4)=sqrt((y0(1)^2)+(y0(2)^2)+(2*(mu1)/r1)+(2*mu/r2)-c);
y0 = [y0(1);y0(2);y0(3); y0(4)]
tspan = [0 20];

options = odeset('RelTol',1e-12,'AbsTol',1e-22,'Events',@events);
             
[t,y,te,ye,ie] = ode45(@f,tspan,y0,options);
q=te(2)
[m,n]=size(y);
xf2=y(m,1);

f2=x1(i)+delx-xf2;
fdash=(f2-f1)/delx;%x1(i)+delx)-xf2)-(x1(i)-xf1), (xf1-xf2)  
e=((x1(i)-xf1)/x1(i));
if(e<0)
    k=e*-1;
else
    k=e;
end
%k=mod(e,1)
if(k<1e-5)
   %uslmanifold(te(2),y0);
   %figure;
   %usrmanifold(te(2),y0)
   %figure;
   %slmanifold(te(2),y0);
   %figure;
   srmanifold(te(2),y0);
   %figure;
   
%plot(y(:,1),y(:,2),'r');
%hillscurve(c);
break;
else
  x1(i+1)= x1(i)-(x1(i)-xf1)/fdash;
  o=x1(i+1);
  end
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

  function [value,isterminal,direction] = events(t,y)
   value = y(2);
    isterminal = 1;        % stop at local minimum
    direction  = 1;         % [local minimum, local maximum]
  end

  % -----------------------------------------------------------------------

end  % 
  
  