function dicompmanofold
% compensated system
ceta=.8;
ts=2;
w= 4/(ceta*ts);
mu = 7.802e-5;
mu1=1-mu;
%load('referencemanifold.mat','y','t')
%xpref=y(:,1);
%ypref=y(:,2);
%xvref=y(:,3);
%yvref=y(:,4);
%xaref=y(:,5);
%yaref=y(:,6);

t0=0;
numSteps=200;
tspan=linspace(0,5,numSteps);
%numerical integration
y0=[0.968371758057117; 0.001299231420374;0.001065939124148;0.016707590741489;0;0;
    0.968360993379390; 0.001300425055020;0.001061146388178;0.016709901591754;0;0;0;0]; %inital condition
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,y0,options);
u1=y(:,13);
u2=y(:,14);
u=sqrt(u1.^2+u2.^2);

figure
plot(t,u1,'r');
grid on

figure
plot(t,u2,'r');
grid on

figure
plot(t,u,'r');
grid on
figure
plot(y(:,1),y(:,2),'r');
hold on;
plot(y(:,7),y(:,8),'b');
p=y
ex=y(1)-y(7);
ey=y(2)-y(8);
vx=y(3)-y(9);
vy=y(4)-y(10);

figure
plot(t(:,1),ex(:,1),'r');
figure
plot(t(:,1),ey(:,1),'b');
figure
plot(t(:,1),vx(:,1),'r');
figure
plot(t(:,1),vy(:,1),'b');
%plot(xpref(:,1),ypref(:,1),'r');
%hold on
%plot(xpsc(:,1),ypsc(:,1),'b');
%controller
%xtr=xasc+2*ux;
%ytr=yasc+2*uy;
%hold on
%plot(xtr(:,1),ytr(:,2),'b');



title('Dynamic Inversion');
ylabel('x,xref');
xlabel('t');



function dydt = f(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
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
      ays
      ux;
      uy];
  end  

end  
  
  