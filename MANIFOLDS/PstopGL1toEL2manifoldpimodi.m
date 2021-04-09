function PstopGL1toEL2manifoldpimodi
% enlarged left unstable manifold of GL1 stopped at poincare section at -pi
mu = 7.802e-5;
mu1=1-mu;
epsilon=2*10^(-5);
haloPeriod= 3.161440780756263;%3.161454133275042;haloPeriod=  2.972432301436252
i0=[ 0.965251234321497;0;0; 0.046053059097282];%[x y xdot ydot]i0=[ 0.969605367097764;0;0;  0.007039341830793];
c= jacobi(i0(1),i0(2),i0(3),i0(4),mu);
%hillscurve(c);
%i0=[ 0.965251237914046;0;0; 0.046053055752129];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
T=3.1*haloPeriod;
k=160;
numSteps=1200;
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
numSteps=1200;
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
tspan=linspace(0, 35,numSteps);
%numerical integration
y0=xU_p(1,:)'; %inital condition
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,y0,options);
hold on;
%plot(y(:,1),y(:,3),'r');
title(c);
xlabel(i0(1));
%plot (mu,0,'*');
%plot (mu1,0,'*');
plot (1.596*mu,0,'*');
plot (1.596*mu1,0,'*');
hold on;
for n=2:k
%compute the initial conditions for the manifold orbits
initialconditions=n;
perturbationVector=monodromy(haloTimes(n),i0)*deltaU;
u_pV=perturbationVector/norm(perturbationVector);
xU_p(n,1:4)=haloPoint(n,1:4)'-epsilon*u_pV;
end
%Compute and plot the unstable manifold
for n=1:k
manifoldLoop=n;
%the positive ones;
t0=0;
numSteps=2000;
tspan=linspace(0,35,numSteps);
%numerical integration
y0=xU_p(n,1:4)'; %inital condition
%options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
%[t,y]=ode113(@f,tspan,y0,options);
options = odeset('RelTol',2.5e-13,'AbsTol',1e-22,'Events',@events1);
             
[t,y,te,ye,ie] = ode113(@f,tspan,y0,options);
plot(y(:,1),y(:,2),'r');
[m,z]=size(y);
yg=y(m,:);
yg1(n,:)=yg;
tg1(n,:)=te;
ye=JEcordinatetransformation(yg,te);
Qe(n,:)=ye;
y=1.596*y;
plot(y(:,1),y(:,2),'r');

end
yg2=Qe
save('ganfile.mat','yg2')
tg2=tg1
save('gantfile.mat','tg2')
%plot(Qe(:,1),Qe(:,3),'r*');
hold on;


[P,Q]=meshgrid(-1.5:0.005:1.5);
r1=((P + mu).^2+Q.^2).^(1/2);
r2=((P + mu-1).^2+Q.^2).^(1/2);
Z=2*[(1/2)*(P.^2 + Q.^2)+(1-mu)./r1+mu./r2];
contour(1.596*P,1.596*Q,1.596*Z,1.596*c,'k')
hold on;

for i=0:-0.005:-2
   plot(i,0,'k-') ;
   hold on;
end




hold on

function dydt = f(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(2))^2);
r2=sqrt((-1+mu+y(1))^2+(y(2))^2);
dydt=[y(3);
      y(4);
      (y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      (y(2)-2*y(3)-(1-mu)*y(2)/(r1)^3-mu*y(2)/(r2)^3)];
  end  

mu =2.523e-5;
mu1=1-mu;
epsilon=2*10^(-5);
haloPeriod= 3.234267818328002;%3.161454133275042;
i0=[1.014011458923152;0;0;0.035771034698967];
c= jacobi(i0(1),i0(2),i0(3),i0(4),mu)%3.00576
hillscurve(c);
%i0=[ 0.965251237914046;0;0; 0.046053055752129];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
T=3.1*haloPeriod;
k=160;
numSteps=1200;
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
numSteps=1200;
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
tspan=linspace(0,50,numSteps);
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
for n=1:k
manifoldLoop=n;
%the positive ones;
t0=0;
numSteps=2000;
tspan=linspace(0,50,numSteps);
%numerical integration
y0=xU_p(n,1:4)'; %inital condition
%options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
%[t,y]=ode113(@f1,tspan,y0,options);
options = odeset('RelTol',2.5e-13,'AbsTol',1e-22,'Events',@events);
             
[t,y,te,ye,ie] = ode113(@f1,tspan,y0,options);
plot(y(:,1),y(:,2),'b');

[m,z]=size(y);
ye=y(m,:);
ye1(n,:)=ye;
te1(n,:)=te;
[m,p]=size(y);
pl1(n,1:4)=y(m,:);

end
ye2=ye1
save('eurfile.mat','ye2')
te2=te1
save('eurtfile.mat','te2')
xlabel('x');
ylabel('xdot');
hold on
L1=pl1;

figure
plot(Qe(:,1),Qe(:,3),'r*');
hold on;
plot(pl1(:,1),pl1(:,3),'b*');
hold on
plot(Qe(:,1),Qe(:,3),'r-');
hold on;
plot(pl1(:,1),pl1(:,3),'b-');
grid on
xlabel('x');
ylabel('xdot');

figure
plot(Qe(:,1),Qe(:,4),'r*');
hold on;
plot(pl1(:,1),pl1(:,4),'b*');
xlabel('x');
ylabel('ydot');
hold on;
plot(Qe(:,1),Qe(:,4),'r-');
hold on;
plot(pl1(:,1),pl1(:,4),'b-');
grid on
xlabel('x');
ylabel('ydot');


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
   value =atan2(y(2),y(1))-1;
    isterminal =1;        % stop at local minimum
    direction  = -1; 
end
function [value,isterminal,direction] = events1(t,y)%stable
   value =atan2(y(2),y(1))-1;
    isterminal =1;        % stop at local minimum
    direction  = -1; 
end
end
