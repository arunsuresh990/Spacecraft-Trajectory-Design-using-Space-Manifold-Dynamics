function unscaledinertialGL1uuslmanifold
%gl1 unstable manifold transformation toinertial unscaled version
%format long
mu = 7.802e-5;
mu1=1-mu;
epsilon=2*10^(-5);
haloPeriod= 2.989309468683970;%3.161454133275042;haloPeriod=  2.972432301436252
i0=[0.968326784284775;0;0;0.016911653732884];%[x y xdot ydot]i0=[ 0.969605367097764;0;0;  0.007039341830793];
c= jacobi(i0(1),i0(2),i0(3),i0(4),mu);
%hillscurve(c);
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
plot (.9706,0,'*');
%plot (1.596*.9706,0,'*');
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
numSteps=200;
tspan=linspace(0,1,numSteps);
%numerical integration
y0=xU_p(1,:)' %inital condition
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,y0,options);
title('GL1 Left Unstable Manifold');
xlabel(mu);
hold on;
for n=2:2
%compute the initial conditions for the manifold orbits
initialconditions=n;
perturbationVector=monodromy(haloTimes(n),i0)*deltaU;
u_pV=perturbationVector/norm(perturbationVector);
xU_p(n,1:4)=haloPoint(n,1:4)'-epsilon*u_pV;
end
%Compute and plot the unstable manifold
for n=2:2
manifoldLoop=n;
%the positive ones;
t0=0;
numSteps=200;
tspan=linspace(0,5,numSteps);
%numerical integration
y0=xU_p(n,1:4)' %inital condition
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,y0,options)

plot(y(:,1),y(:,2),'k');
[m,p]=size(y);
[m,z]=size(y);
fg=y(m,:);
p=y
rg=1.070e6;
vg=10.909;
wg=1.02e-5;
Tg=7.15*24*60*60;

pr(:,1)=y(:,1)*rg;
pr(:,2)=y(:,2)*rg;
%pr(:,3)=0*rg;
vr(:,1)=y(:,3)*vg;
vr(:,2)=y(:,4)*vg;
%vr(:,3)=0*vg;
tr=(Tg/2*pi)*t;

Pin=pr;
vg=vr;
tg=tr;
end


rg=1.070e6;
wg=1.02e-5;
re=6.711e5;
we=2.047e-5;
Tg=7.15*24*60*60;
Te=3.55*24*60*60;




pos=Pin;
%yi=Pi(:,2)


figure
%ganymede to inertial
plot(0,0,'k*');
hold on;
plot(Pin(:,1),Pin(:,2),'r');
t=tg;
for i=1:numSteps
theetag(i,1)=(wg*t(i));
TG2I=[cos(theetag(i)) -sin(theetag(i));sin(theetag(i)) cos(theetag(i))];
Pgi(i,:)=TG2I*(Pin(i,:))';
end

o1=Pgi;

figure
plot(Pgi(:,1),Pgi(:,2),'b');
hold on;

%inertial ganymede to Ganymede rotating frame
for i=1:numSteps
theetag(i,1)=(wg*t(i));
c=Pgi(i,:)';
TI2G=[cos(theetag(i)) sin(theetag(i));-sin(theetag(i)) cos(theetag(i))];
Pgg(i,:)=TI2G*c;
end
o3=Pgg;

%plot(Pig(:,1),Pig(:,2),'r');


%inertial Ganymede to Europa rotating frame
%numSteps=200;
%t=linspace(0,Tg/4,numSteps)';
for i=1:numSteps
theetag(i,1)=(we*t(i));
%TG2I=[cos(theetag(i)) -sin(theetag(i));sin(theetag(i)) cos(theetag(i))];
%Pgi(i,:)=TG2I*Pg;
c=Pgi(i,:)';
TI2E=[cos(theetag(i)) sin(theetag(i));-sin(theetag(i)) cos(theetag(i))];
Pge(i,:)=TI2E*c;
end
o4=Pge;





figure
title('Ganymede Rotating Frame');
plot(0,0,'k*');
hold on;
plot(Pgg(:,1),Pgg(:,2),'r');


figure
title('Europa Rotating Frame');
plot(0,0,'k*');
hold on;
plot(Pge(:,1),Pge(:,2),'b');



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