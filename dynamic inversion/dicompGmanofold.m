function dicompGmanofold
% compensated system
ceta=.8;
ts=2;
w= 4/(ceta*ts);
mu = 7.802e-5;
mu1=1-mu;

t0=0;
numSteps=1000;
dur=5;
tspan=linspace(0,dur,numSteps);
%%numerical integration, first 6 entries is position velocity and
%acceleration of REFERENCE manifold
y0=[0.967071758057117; 0.001299231420374;0.001065939124148;0.016707590741489;0;0;
    0.968360993379390; 0.001300425055020;0.001061146388178;0.016709901591754;0;0]; %inital condition
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,y0,options);
plot(y(:,1),y(:,2),'r');
hold on;
plot(y(:,7),y(:,8),'b');
hleg1 = legend('Reference','Actual');
title('Compensated System after Dynamic Inversion ');
ylabel('x,xref');
xlabel('t');

%ERRORS
ex(:,1)=y(:,1)-y(:,7);
ey(:,1)=y(:,2)-y(:,8);
vx(:,1)=y(:,3)-y(:,9);
vy(:,1)=y(:,4)-y(:,10);

% calling uncompensated system
[y1,t1]=diuncompGmanofold(y0,numSteps,dur);
ex1(:,1)=y1(:,1)-y1(:,7);
ey1(:,1)=y1(:,2)-y1(:,8);
vx1(:,1)=y1(:,3)-y1(:,9);
vy1(:,1)=y1(:,4)-y1(:,10);

%Error plots
figure
plot(t,ex,'r-');
hold on
plot(t1,ex1,'b-');
title('Error in position in x direction');
ylabel('Error in x');
xlabel('t');

figure
plot(t,ey,'r-');
hold on
plot(t1,ey1,'b-');
title('Error in position in y direction');
ylabel('Error in y');
xlabel('t');

figure
plot(t,vx,'r-');
hold on
plot(t1,vx1,'b-');
title('Error in velocity in x direction');
ylabel('Error in x velocity');
xlabel('t');

figure
plot(t,vy,'r-');
hold on
plot(t1,vy1,'b-');
title('Error in velocity in y direction');
ylabel('Error in y velocity');
xlabel('t');



function dydt = f(t,y)
r1r=sqrt((mu+y(1))^2+(y(2))^2);
r2r=sqrt((-1+mu+y(1))^2+(y(2))^2);
r1s=sqrt((mu+y(7))^2+(y(8))^2);
r2s=sqrt((-1+mu+y(7))^2+(y(8))^2);
axr=y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1r)^3-mu/r2r^3*(y(1)-(1-mu));
ayr=y(2)-2*y(3)-(1-mu)*y(2)/(r1r)^3-mu*y(2)/(r2r)^3;
ux=(1/2)*((axr-(y(7)+2*y(10)-(1-mu)*(y(1)+mu)/(r1s)^3-mu/r2s^3*(y(7)-(1-mu))))-2*ceta*w*(y(9)-y(3))-w^2*(y(7)-y(1)));
uy=(1/2)*((ayr-(y(8)-2*y(9)-(1-mu)*y(8)/(r1s)^3-mu*y(8)/(r2s)^3))-2*ceta*w*(y(10)-y(4))-w^2*(y(8)-y(2)));
axs=y(7)+2*y(10)-(1-mu)*(y(7)+mu)/(r1s)^3-mu/r2s^3*(y(7)-(1-mu))+2*ux;
ays=y(8)-2*y(9)-(1-mu)*y(8)/(r1s)^3-mu*y(8)/(r2s)^3+2*uy;
dydt=[y(3);
      y(4);
      (y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1r)^3-mu/r2r^3*(y(1)-(1-mu)));
      (y(2)-2*y(3)-(1-mu)*y(2)/(r1r)^3-mu*y(2)/(r2r)^3);
      axr;
      ayr;
      y(9);
      y(10);
      (y(7)+2*y(10)-(1-mu)*(y(1)+mu)/(r1s)^3-mu/r2s^3*(y(7)-(1-mu))+2*ux);
      (y(8)-2*y(9)-(1-mu)*y(8)/(r1s)^3-mu*y(8)/(r2s)^3+2*uy);
      axs;
      ays];
  end  

end  
  
  