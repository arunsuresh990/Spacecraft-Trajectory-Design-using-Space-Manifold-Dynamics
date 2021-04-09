function simpletrajectoryode
% Problem parameters
mu = 7.802e-5;
mu1=1-mu;
c= 3.0058;
%hillscurve(c);
y0(1)=   0.978332646874192;
y0(2)= -0.008036026547563  ;
y0(3)=0.009858514255237  ;
%y0=[    0.968336849830872 ; 0 ; 0;  0.016911577500865];
r1=sqrt((mu+y0(1))^2+(y0(3))^2);
r2=sqrt((-mu1+y0(1))^2+(y0(3))^2);
y0(4)=sqrt((y0(1)^2)+(y0(3)^2)+(2*(mu1)/r1)+(2*mu/r2)-c)
%numerical integration split up into a number of  steps
y0 = [y0(1);y0(2);y0(3); y0(4)];
numSteps=200;
p=1
ti=0;
tend=1;
k=tend/numSteps

for i=0:k:tend-k
x(p,:)=y0;   
tf=ti+k;
tspan = [ti tf];
options = odeset('RelTol',1e-12,'AbsTol',1e-8);             
[t,y] = ode45(@f,tspan,y0,options);
[m,z]=size(y);
y0=y(m,:);
p=p+1;
plot(y(:,1),y(:,2),'r');
hold on
plot (mu1,0,'*');
ti=tf;
end
l=x
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

  
  