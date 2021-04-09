function transformationmanifold2
%transform  Gl1 unstable manifold in JG rotating frame to JE rotating frame

mu = 7.802e-5;
mu1=1-mu;
epsilon=2*10^(-5);
haloPeriod= 3.161440780756263;%3.161454133275042;haloPeriod=  2.972432301436252
i0=[ 0.965251234321497;0;0; 0.046053059097282];%[x y xdot ydot]i0=[ 0.969605367097764;0;0;  0.007039341830793];
c= jacobi(i0(1),i0(2),i0(3),i0(4),mu);
hillscurve(c);
hold on
C=3.0058;
scaledhillscurve(C);



%i0=[ 0.965251237914046;0;0; 0.046053055752129];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
T=3.1*haloPeriod;
k=40;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);
%t=haloPeriod;
%tspan = [0 t];
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,i0,options);%FUNCTION CALL
%plot(y(:,1),y(:,2),'k.');
hold on;
%plot (.9706,0,'*');
plot (1.596*.9706,0,'*');
mon=monodromy(haloPeriod,i0);%FUNCTION CALL
o=det(mon);
[V,E] = eigs(mon);
moneig=diag(E);
%get the parameterization points of the periodic orbit
numSteps=300;
h=round((numSteps-1)/k);
for n=1:k
haloPoint(n,1:4)=y((n-1)*h+1,1:4);
haloTimes(n)=t((n-1)*h+1);
end


I=eye(4);
%The un-stable eigenvector is
deltaU=V(1:4,1)
if (deltaU(1)<0)
    deltaU=-deltaU
end

%Compute an orbit on the unstable manifold
%begining at the origin of the orbit
xU_p(1,1:4)=haloPoint(1,1:4)'-epsilon*I*deltaU;
%compute the unstable orbit
%the positive ones;
t0=0;
numSteps=2000;
tspan=linspace(0, 25,numSteps);
%numerical integration
y0=xU_p(1,:)'; %inital condition
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,y0,options);
hold on;
%plot(y(:,1),y(:,3),'r');
title(c);
xlabel(i0(1));
plot (1.596*mu,0,'*');
plot (1.596*mu1,0,'*');
%plot (mu,0,'*');
%plot (mu1,0,'*');
hold on;
for n=2:k
%compute the initial conditions for the manifold orbits
initialconditions=n;
perturbationVector=monodromy(haloTimes(n),i0)*deltaU;
u_pV=perturbationVector/norm(perturbationVector);
xU_p(n,1:4)=haloPoint(n,1:4)'-epsilon*u_pV;
end
%Compute and plot the unstable manifold
for n=1:1
manifoldLoop=n;
%the positive ones;
t0=0;
numSteps=200;
tspan=linspace(0,25,numSteps);
%numerical integration
y0=xU_p(n,1:4)'; %inital condition
%options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
%[t,y]=ode113(@f,tspan,y0,options);
hold on;
options = odeset('RelTol',1e-12,'AbsTol',1e-22,'Events',@events);
             
[t,y,te,ye,ie] = ode113(@f,tspan,y0,options);
[m,z]=size(y);
for i=1:m
yg=y(i,:);
tg=t(i);
tg1(i,:)=tg;
yg1(i,:)=yg;
ye=JEcordinatetransformation(yg,tg);
Qe(i,:)=ye;
end
y1=yg1
t1=tg1
Q1=Qe
l=Q1(:,1);
xdot(1)=0;
for i=2:m
xdot(i,1) =(l(i)-l(i-1))/(t(i)-t(i-1));
end
h=xdot
l=Q1(:,2);
ydot(1)=0;
for i=2:m
ydot(i,1) =(l(i)-l(i-1))/(t(i)-t(i-1));
end
j=ydot
plot(y(:,1),y(:,2),'b');
figure
t = linspace(0,2*pi);
plot(1.596*cos(t),1.596*sin(t),'k*');
hold on
plot(Qe(:,1),Qe(:,2),'r');
hold on;
end
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
function [value,isterminal,direction] = events(t,y)
   
    
   value =atan2(y(2),y(1))-1;
    isterminal =1;        % stop at local minimum
    direction  = -1;         % [local minimum, local maximum]
end
end