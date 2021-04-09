function PstopGL2toGL1manifold1mu
%left unstable of L2 and right stable of L1 trajectories cut at poincare
%section at mu
format long;
mu = 7.802e-5;
mu1=1-mu;
epsilon=2*10^(-4);
haloPeriod=3.181963172624672;
i0=[ 1.022851960426745;0;0; 0.039327988046747];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
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


xU_p(1,1:4)=haloPoint(1,1:4)'-epsilon*I*deltaU;

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
xU_p(n,1:4)=haloPoint(n,1:4)'-epsilon*u_pV;
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


haloPeriod= 3.075183533429036;
i0=[ 0.966219738362142;0;0; 0.036117254575734];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
c1= jacobi(i0(1),i0(2),i0(3),i0(4),mu)%3.00576
hillscurve(c1);
T=3.1*haloPeriod;
k=40;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);
%t=haloPeriod;
%tspan = [0 t];
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t1,y1]=ode113(@f,tspan,i0,options);%FUNCTION CALL
plot(y1(:,1),y1(:,2),'k.');
hold on;
plot (.9706,0,'*');
mon=monodromy(haloPeriod,i0);%FUNCTION CALL
[V,E] = eigs(mon);
moneig=diag(E);
%get the parameterization points of the periodic orbit
numSteps=300;
h=round((numSteps-1)/k);
for n=1:k
haloPoint(n,1:4)=y1((n-1)*h+1,1:4);
haloTimes(n)=t1((n-1)*h+1);
end


I=eye(4);
%The stable eigenvector is
deltaS=V(1:4,4);
if (deltaS(1)<0)
   deltaS=-deltaS;
end
%Compute an orbit on the unstable manifold
%begining at the origin of the orbit
xU_p(1,1:4)=haloPoint(1,1:4)'+epsilon*I*deltaS;
%compute the unstable orbit
%the positive ones;
t0=0;
numSteps=2000;
tspan=linspace(0,5,numSteps);
%numerical integration
y0=xU_p(1,:)'; %inital condition
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t1,y1]=ode113(@f1,tspan,y0,options);
%plot(y(:,1),y(:,3),'r');
plot (mu,0,'*');
plot (mu1,0,'*');
hold on;
for n=2:k
%compute the initial conditions for the manifold orbits
initialconditions=n;
perturbationVector=monodromy(haloTimes(n),i0)*deltaS;
u_pV=perturbationVector/norm(perturbationVector);
xU_p(n,1:4)=haloPoint(n,1:4)'+epsilon*u_pV;
end
%Compute and plot the unstable manifold
for n=2:k
manifoldLoop=n;
%the positive ones;
t0=0;
numSteps=2000;
tspan=linspace(0, 5,numSteps);
%numerical integration
y0=xU_p(n,1:4)'; %inital condition
options = odeset('RelTol',2.5e-13,'AbsTol',1e-22,'Events',@events1);
             
[t1,y1,te1,ye1,ie1] = ode113(@f1,tspan,y0,options);

%options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
%[t,y]=ode113(@f1,tspan,y0,options);
%plot the unstable pos orbits as red lines
plot(y1(:,1),y1(:,2),'b');
hold on;
grid on;
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
function [value,isterminal,direction] = events(t,y)%stable
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

