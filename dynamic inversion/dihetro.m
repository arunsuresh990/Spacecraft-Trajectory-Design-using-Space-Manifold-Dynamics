function dihetro
% compensated system
ceta=.8;
ts=2;
w= 4/(ceta*ts);
mu = 7.802e-5;
mu1=1-mu;

t0=0;
numSteps=1000;
dur=3;
tspan=linspace(0,dur,numSteps);
%%numerical integration, first 6 entries is position velocity and
%acceleration of REFERENCE manifold
y0=[1.034731870033446 ;  0.028798141703113;  -0.009107415607789 ; -0.013919377706517;0;0;
 1.039731870033446 ;  0.028398141703113;  -0.009107415607789 ; -0.013919377706517;0;0;0;0];
 %inital condition %inital condition 6 for reference trajectory, 6 for actual, 2 for control vector
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t1,y1]=ode113(@f,tspan,y0,options)
[m,z]=size(t1);
t01=t1(m,:);
[y1u,t1u]=diuncompEmanofold(y0,numSteps,dur);
[m,z]=size(t1u);
t01u=t1u(m,:);
t0=0;
numSteps=1000;
dur=3;
tspan=linspace(0,dur,numSteps);
%%numerical integration, first 6 entries is position velocity and
%acceleration of REFERENCE manifold
y0=[0.99992198;.01557542;- 0.056623322622419;-.00037135;0;0;
 0.99992198;.01557542;- 0.069623322622419;-.00197135;0;0;0;0];
 %inital condition %inital condition 6 for reference trajectory, 6 for actual, 2 for control vector
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t2,y2]=ode113(@f,tspan,y0,options)

t2=t2+t01;
y=[y1;y2]
t=[t1;t2]

[y2u,t2u]=diuncompEmanofold(y0,numSteps,dur);
t2u=t2u+t01u;
y1=[y1u;y2u]
t1=[t1u;t2u]

u1=y(:,13);
u2=y(:,14);
u=sqrt(u1.^2+u2.^2);

plot(t,u,'r');
title('Control Signal ');
ylabel('Control Signal ');
xlabel('t');
grid on

i=1;
if i==1
    
figure
plot(t,u1,'k');
title('Control Signal in x direction');
ylabel('Control Signal in x');
xlabel('t');
grid on

figure
plot(t,u2,'k');
title('Control Signal in y direction');
ylabel('Control Signal in y');
xlabel('t');
grid on

figure
plot(y(:,1),y(:,2),'k');
hold on;
plot(y(:,7),y(:,8),'k--');
hleg1 = legend('Reference','Actual');
title('Compensated System after Dynamic Inversion ');
ylabel('x,xref');
xlabel('t');

%ERRORS
ex(:,1)=y(:,1)-y(:,7);
ey(:,1)=y(:,2)-y(:,8);
vx(:,1)=y(:,3)-y(:,9);
vy(:,1)=y(:,4)-y(:,10);
ax(:,1)=y(:,5)-y(:,11);
ay(:,1)=y(:,6)-y(:,12);

% calling uncompensated system
%[y1,t1]=diuncompEmanofold(y0,numSteps,dur);
ex1(:,1)=y1(:,1)-y1(:,7);
ey1(:,1)=y1(:,2)-y1(:,8);
vx1(:,1)=y1(:,3)-y1(:,9);
vy1(:,1)=y1(:,4)-y1(:,10);
ax1(:,1)=y1(:,5)-y1(:,11);
ay1(:,1)=y1(:,6)-y1(:,12);

%Error plots
figure
plot(t,ex,'k');
hold on
plot(t1,ex1,'k--');
hleg1 = legend('With DI','Without DI');
title('Error in position in x direction');
ylabel('Error in x');
xlabel('t');
grid on;

figure
plot(t,ey,'k');
hold on
plot(t1,ey1,'k--');
hleg1 = legend('With DI','Without DI');
title('Error in position in y direction');
ylabel('Error in y');
xlabel('t');
grid on;

figure
plot(t,vx,'k');
hold on
plot(t1,vx1,'k--');
hleg1 = legend('With DI','Without DI');
title('Error in velocity in x direction');
ylabel('Error in x velocity');
xlabel('t');
grid on;

figure
plot(t,vy,'k');
hold on
plot(t1,vy1,'k--');
hleg1 = legend('With DI','Without DI');
title('Error in velocity in y direction');
ylabel('Error in y velocity');
xlabel('t');
grid on;

figure
plot(t,ax,'r-');
hold on
plot(t1,ax1,'b-');
title('Error in acceleration in x direction');
ylabel('Error in x acceleration');
xlabel('t');
grid on;

figure
plot(t,ay,'r-');
hold on
plot(t1,ay1,'b-');
title('Error in acceleration in y direction');
ylabel('Error in y acceleration');
xlabel('t');
grid on;
end


function dydt = f(t,y)
   
r1r=sqrt((mu+y(1))^2+(y(2))^2);
r2r=sqrt((-1+mu+y(1))^2+(y(2))^2);
r1s=sqrt((mu+y(7))^2+(y(8))^2);
r2s=sqrt((-1+mu+y(7))^2+(y(8))^2);
axr=y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1r)^3-mu/r2r^3*(y(1)-(1-mu));
ayr=y(2)-2*y(3)-(1-mu)*y(2)/(r1r)^3-mu*y(2)/(r2r)^3;
ux=(1/2)*((axr-(y(7)+2*y(10)-(1-mu)*(y(7)+mu)/(r1s)^3-mu/r2s^3*(y(7)-(1-mu))))-2*ceta*w*(y(9)-y(3))-w^2*(y(7)-y(1)));
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
      ays;
      ux;
      uy];
    

end
end