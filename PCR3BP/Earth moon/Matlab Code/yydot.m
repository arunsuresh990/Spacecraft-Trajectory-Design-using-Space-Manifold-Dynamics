function yydot
%right stable manifold of moon L2 stopped at poincare section at -mu(earth
%position)
mu=1.215e-2;
mu1=1-mu;
epsilon=2*10^(-4);
haloPeriod= 3.397804099470672;%3.161454133275042;haloPeriod=  2.972432301436252
i0=[1.129091528063946;0;0;0.134455226022298];%[x y xdot ydot]i0=[ 0.969605367097764;0;0;  0.007039341830793];
c= jacobi(i0(1),i0(2),i0(3),i0(4));
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
%plot(y(:,1),y(:,2),'k.');
hold on;
mon=monodromy(haloPeriod,i0)%FUNCTION CALL
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
numSteps=2000;
tspan=linspace(0,30,numSteps);
%numerical integration
y0=xU_p(1,:)'; %inital condition
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f1,tspan,y0,options);
%plot(y(:,1),y(:,3),'r');
title(c);
xlabel(i0(1));
hold on;
for n=2:k
%compute the initial conditions for the manifold orbits
initialconditions=n
perturbationVector=monodromy(haloTimes(n),i0)*deltaS;
u_pV=perturbationVector/norm(perturbationVector);
xU_p(n,1:4)=haloPoint(n,1:4)'+epsilon*u_pV;
end
%Compute and plot the unstable manifold
for n=1:k
manifoldLoop=n
%the positive ones;
t0=0;
numSteps=2000;
tspan=linspace(0,30,numSteps);
%numerical integration
y0=xU_p(n,1:4)'; %inital condition
options = odeset('RelTol',2.5e-13,'AbsTol',1e-22,'Events',@events);
             
[t,y,te,ye,ie] = ode113(@f1,tspan,y0,options);
[m,z]=size(y);
plot(y(:,1),y(:,2),'b');
hold on;
yg=y(m,:);
yg1(n,:)=yg;
tg1(n,:)=te;
ye=MEcordinatetransformation(yg,te);
Qe(n,:)=ye;
end
figure
plot(Qe(:,2),Qe(:,4),'r-');
hold on
load('yydotSE.mat','yg1') ;
plot(yg1(:,2),yg1(:,4),'b-');





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
   value =atan2(y(2),y(1))-2.25;
    isterminal =1;        % stop at local minimum
    direction  = 0;         % [local minimum, local maximum]
 end

end