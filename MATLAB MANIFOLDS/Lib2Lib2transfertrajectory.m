function L2L2transfertrajectory%left unstable of L2 and right stable of L1
format long;
mu = 7.802e-5;
mu1=1-mu;
epsilon=2*10^(-5);
haloPeriod=3.106407530366246;
i0=[ 1.027613550302901;0;0;0.013486361810618];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
c= jacobi(i0(1),i0(2),i0(3),i0(4),mu)%3.00576
hillscurve(c);
k=40;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,i0,options);%FUNCTION CALL
plot(y(:,1),y(:,2),'k.');
hold on;
plot (1.0298,0,'*');
mon=monodromy(haloPeriod,i0);%FUNCTION CALL
[V,E] = eigs(mon);
moneig=diag(E);
numSteps=300;
h=round((numSteps-1)/k);
for n=1:k
haloPoint(n,1:4)=y((n-1)*h+1,1:4);
haloTimes(n)=t((n-1)*h+1);
end


I=eye(4);
deltaU=V(1:4,1);
if (deltaU(1)<0)
    deltaU=-deltaU;
end


xU_p(1,1:4)=haloPoint(1,1:4)'+epsilon*I*deltaU;

t0=0;
numSteps=2000;
tspan=linspace(0, 5,numSteps);
y0=xU_p(1,:)'; 
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,y0,options);
title(c);
xlabel(i0(1));
plot (mu,0,'*');
plot (mu1,0,'*');
hold on;
for n=2:k
initialconditions=n;
perturbationVector=monodromy(haloTimes(n),i0)*deltaU;
u_pV=perturbationVector/norm(perturbationVector);
xU_p(n,1:4)=haloPoint(n,1:4)'+epsilon*u_pV;
end
for n=2:k
manifoldLoop=n;
t0=0;
numSteps=2000;
tspan=linspace(0, 5,numSteps);
y0=xU_p(n,1:4)'; %inital condition
%options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
%[t,y]=ode113(@f,tspan,y0,options);
options = odeset('RelTol',2.5e-13,'AbsTol',1e-22,'Events',@events);
             
[t,y,te,ye,ie] = ode113(@f,tspan,y0,options);

plot(y(:,1),y(:,2),'r');
axis([.96 1.052 -.04 .04]);
for i=-0.02:0.0005:0.02
   plot(mu1,i,'k.') ;
   hold on;
end
end


haloPeriod=3.106407530366246;
i0=[ 1.027613550302901;0;0;0.013486361810618];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
c= jacobi(i0(1),i0(2),i0(3),i0(4))%3.00576
hillscurve(c);
k=40;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,i0,options);%FUNCTION CALL
plot(y(:,1),y(:,2),'k.');
hold on;
plot (1.0298,0,'*');
mon=monodromy(haloPeriod,i0);%FUNCTION CALL
[V,E] = eigs(mon);
moneig=diag(E);
numSteps=300;
h=round((numSteps-1)/k);
for n=1:k
haloPoint(n,1:4)=y((n-1)*h+1,1:4);
haloTimes(n)=t((n-1)*h+1);
end


I=eye(4);
deltaU=V(1:4,4);
if (deltaU(1)<0)
    deltaU=-deltaU;
end


xU_p(1,1:4)=haloPoint(1,1:4)'+epsilon*I*deltaU;

t0=0;
numSteps=2000;
tspan=linspace(0, 5,numSteps);
y0=xU_p(1,:)'; 
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,y0,options);
title(c);
xlabel(i0(1));
plot (mu,0,'*');
plot (mu1,0,'*');
hold on;
for n=2:k
initialconditions=n;
perturbationVector=monodromy(haloTimes(n),i0)*deltaU;
u_pV=perturbationVector/norm(perturbationVector);
xU_p(n,1:4)=haloPoint(n,1:4)'+epsilon*u_pV;
end
for n=2:k
manifoldLoop=n;
t0=0;
numSteps=2000;
tspan=linspace(0, 5,numSteps);
y0=xU_p(n,1:4)'; %inital condition
%options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
%[t,y]=ode113(@f,tspan,y0,options);
options = odeset('RelTol',2.5e-13,'AbsTol',1e-22,'Events',@events);
             
[t,y,te,ye,ie] = ode113(@f,tspan,y0,options);

plot(y(:,1),y(:,2),'r');
axis([.96 1.052 -.04 .04]);

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
   value = y(1)-(1-mu);
    isterminal =1;        % stop at local minimum
    direction  = -1;         % [local minimum, local maximum]
end
 function [value,isterminal,direction] = events1(t1,y1)
   value = y1(1)-(1-mu);
    isterminal =1;        % stop at local minimum
    direction  = 1;         % [local minimum, local maximum]
 end
end
