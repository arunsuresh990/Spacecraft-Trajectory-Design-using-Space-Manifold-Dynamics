function EL2ssrmanifold
% right stable manifold of EL2 stopped at poincare section at -pi
mu =2.523e-5;
mu1=1-mu;
epsilon=2*10^(-5);
haloPeriod= 3.234267818328002;%3.161454133275042;
i0=[1.014011458923152;0;0;0.035771034698967];
c= jacobi(i0(1),i0(2),i0(3),i0(4),mu)%3.00576
hillscurve(c);
%i0=[ 0.965251237914046;0;0; 0.046053055752129];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
T=3.1*haloPeriod;
k=40;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);
%t=haloPeriod;
%tspan = [0 t];
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,i0,options);%FUNCTION CALL
plot(y(:,1),y(:,2),'k.');
hold on;
plot (1.020,0,'*');
mon=monodromyE(haloPeriod,i0)%FUNCTION CALL
o=det(mon)
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
%The stable eigenvector is
deltaS=V(1:4,4)
if (deltaS(1)<0)
   deltaS=-deltaS
end
%Compute an orbit on the unstable manifold
%begining at the origin of the orbit
xU_p(1,1:4)=haloPoint(1,1:4)'+epsilon*I*deltaS;
%compute the unstable orbit
%the positive ones;
t0=0;
numSteps=200;
tspan=linspace(0,20,numSteps);
%numerical integration
y0=xU_p(1,:)'; %inital condition
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f1,tspan,y0,options);
%plot(y(:,1),y(:,3),'r');
title(c);
xlabel(i0(1));
plot (-mu,0,'*');
plot (mu1,0,'*');
hold on;
for n=2:k
%compute the initial conditions for the manifold orbits
initialconditions=n;
perturbationVector=monodromyE(haloTimes(n),i0)*deltaS;
u_pV=perturbationVector/norm(perturbationVector);
xU_p(n,1:4)=haloPoint(n,1:4)'+epsilon*u_pV;
end
%Compute and plot the unstable manifold
for n=2:2
manifoldLoop=n;
%the positive ones;
t0=0;
numSteps=2000;
tspan=linspace(0, 50,numSteps);
%numerical integration
y0=xU_p(n,1:4)'; %inital condition
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f1,tspan,y0,options);
plot(y(:,1),y(:,2),'b');
[m,z]=size(y);
fg=y(m,:)
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
end