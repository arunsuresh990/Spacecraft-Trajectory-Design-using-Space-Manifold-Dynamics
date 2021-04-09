function GanymedetoEuropa
format long;
rg=1.070e6;
re=6.711e5;
Tg=7.15*24*60*60;
mu = 7.802e-5;
mu1=1-mu;
i0=[0.615742937358920 ;  0.468137770650280 ; -0.247571270026253  ; 0.340976378442258];%obtained from poincare section and transformed back to ganymede frame
numSteps=2000;
tspan=linspace(0,5,numSteps);
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f1,tspan,i0,options);
[yi,ye]=GframetoEframe(y,t,numSteps);
tr=(Tg/2*pi)*t;
plot(tr(:,1),yi(:,1),'b');
figure
plot(yi(:,1),yi(:,2),'r');

hold on;


%ganymedetrajectory;
hold on;
t = linspace(0,2*pi);
plot(re*cos(t),re*sin(t),'k*');
hold on
t = linspace(0,2*pi);
plot(rg*cos(t),rg*sin(t),'k*');
M=(1-mu);
J=-mu;
plot(J,0,'r*');
hold on;
plot(re*M,0,'r*');
hold on;
plot(rg*M,0,'r*');



mu =2.523e-5;
mu1=1-mu;
epsilon=2*10^(-5);
haloPeriod= 3.234267818328002;%3.161454133275042;
i0=[0.978970256289329  ; 0.669749586599251 ;  0.177315345423028 ; -0.24256920193660];%0.27467,0.361687574
%i01=[-1.222124794736756  -0.000000000000059  -0.003918596097935   0.355638203395082];
c= 3.0028
%hillscurve(c);
numSteps=2000;
tspan=linspace(0,5,numSteps);
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y2]=ode113(@f,tspan,i0,options);

[yi,yg]=EframetoGframe(y2,t,numSteps);

plot(yi(:,1),yi(:,2),'b');
hold on;
figure
plot(y2(:,1),y2(:,2),'b');
hold on;
plot(ye(:,1),ye(:,2),'r');

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
end