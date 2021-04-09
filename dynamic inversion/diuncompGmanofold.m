function [y1,t1]= diuncompGmanofold(i0,numsteps,dur)
% uncompensated system
ceta=.5;
ts=10;
w= 10/(ceta*ts);
mu = 7.802e-5;
mu1=1-mu;

t0=0;
numSteps=numsteps;
tspan=linspace(0,dur,numSteps);
%numerical integration
%y0=[0.968371758057117; 0.001299231420374;0.001065939124148;0.016707590741489;0;0;
   % 0.968360993379390; 0.001300425055020;0.001061146388178;0.016709901591754;0;0]; 
y0=i0;%inital condition
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,y0,options);
figure
plot(y(:,1),y(:,2),'r');
hold on;
plot(y(:,7),y(:,8),'b');
hleg1 = legend('Reference','Actual');
y1=y;
t1=t;
ex(:,1)=y(:,1)-y(:,7);
ey(:,1)=y(:,2)-y(:,8);
vx(:,1)=y(:,3)-y(:,9);
vy(:,1)=y(:,4)-y(:,10);



%figure
%grid on
%plot(t,ex,'r-');


%plot(xpref(:,1),ypref(:,1),'r');
%hold on
%plot(xpsc(:,1),ypsc(:,1),'b');
%controller
%xtr=xasc+2*ux;
%ytr=yasc+2*uy;
%hold on
%plot(xtr(:,1),ytr(:,2),'b');

%save('uncompmanifold.mat','y','t');

title('Uncompensated System before Dynamic Inversion');
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
axs=y(7)+2*y(10)-(1-mu)*(y(7)+mu)/(r1s)^3-mu/r2s^3*(y(7)-(1-mu));
ays=y(8)-2*y(9)-(1-mu)*y(8)/(r1s)^3-mu*y(8)/(r2s)^3;
%ux=(1/2)*((axr-(y(7)+2*y(10)-(1-mu)*(y(1)+mu)/(r1s)^3-mu/r2s^3*(y(7)-(1-mu))))-2*ceta*w*(y(9)-y(3))-w^2*(y(7)-y(1)));
%uy=(1/2)*((axr-(y(8)-2*y(9)-(1-mu)*y(8)/(r1s)^3-mu*y(8)/(r2s)^3))-2*ceta*w*(y(9)-y(3))-w^2*(y(7)-y(1)));
dydt=[y(3);
      y(4);
      (y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1r)^3-mu/r2r^3*(y(1)-(1-mu)));
      (y(2)-2*y(3)-(1-mu)*y(2)/(r1r)^3-mu*y(2)/(r2r)^3);
      axr;
      ayr;
      y(9);
      y(10);
      (y(7)+2*y(10)-(1-mu)*(y(1)+mu)/(r1s)^3-mu/r2s^3*(y(7)-(1-mu)));
      (y(8)-2*y(9)-(1-mu)*y(8)/(r1s)^3-mu*y(8)/(r2s)^3);
      axs;
      ays];
  end  

end  
  
  