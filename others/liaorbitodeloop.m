function liaorbitodeloop
% Problem parameters
mu = 7.802e-5;
mu1=1-mu;
c= 3.0058;%3.00576%3.0075940671432188;
L=0.97058;
p1=L-.003;
for i=1:3
y0(1)=p1 ;%0.96599912805843%0.969539872159348, 1.0201059805843
y0(3)=0;
y0(2)=0;
%y0(4)=0.016911657905671;
r1=sqrt((mu+y0(1))^2+(y0(3))^2);
r2=sqrt((-mu1+y0(1))^2+(y0(3))^2);
y0(4)=sqrt((y0(1)^2)+(y0(3)^2)+(2*(mu1)/r1)+(2*mu/r2)-c);%0.007530309627883
y0 = [y0(1);y0(2);y0(3); y0(4)]
tspan = [0 40];

options = odeset('RelTol',1e-12,'AbsTol',1e-8,'Events',@events1);             
[t,y,te,ye,ie] = ode45(@f,tspan,y0,options);
[m1,n1]=size(y);
x1=y(m1,1);
figure
  plot(y(:,1),y(:,3),'r');

options = odeset('RelTol',1e-12,'AbsTol',1e-8,'Events',@events2);             
[t,x,te,ye,ie] = ode45(@f,tspan,y0,options);
[m2,n2]=size(y);
x2=y(m2,1);
 figure
 plot(x(:,1),x(:,3),'b');
k1=x1-x2
if(k1>0)
    C1=c;
    init1=p1;
    plot(y(:,1),y(:,3),'r');
    hold on;
    plot (0.97058,0,'*');
    %liaorbitwwithnrloop(C1,init1);
    break
    
else
    k2=x1-L
     b=(p1+L)/2;
     a=p1-L;
    if(k2>0) 
       p2=b-abs(a);
   
    else
       p2=b;
       
    end     
    end
p1=p2;       
end     
  function dydt = f(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(3))^2);
r2=sqrt((-1+mu+y(1))^2+(y(3))^2);
dydt=[y(2);
     (y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      y(4);
     (y(3)-2*y(2)-(1-mu)*y(3)/(r1)^3-mu*y(3)/(r2)^3)];
  end  

  % -----------------------------------------------------------------------

  function [value,isterminal,direction] = events1(t,y)
   value = y(3);
    isterminal = 1;        % stop at local minimum
    direction  = -1;         
  end

function [value,isterminal,direction] = events2(t,y)
   value = y(3);
    isterminal = 1;        % stop at local minimum
    direction  = 1;         
  end

  % -----------------------------------------------------------------------

end  % liaorbitode

  
  