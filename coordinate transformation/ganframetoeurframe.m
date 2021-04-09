function [yi,ye]= ganframetoeurframe(yg,tg,num)
%function ganframetoeurframe
%mu = 7.802e-5;
%mu1=1-mu;
%epsilon=2*10^(-5);
numSteps=num;
tspan=linspace(0,13,numSteps);
%y0=[-0.768579953607849;-0.000000000000006;-0.007378129569614;-0.432726026515497];%[x y xdot ydot]i0=[ 0.969605367097764;0;0;  0.007039341830793];%inital condition
%options=odeset('RelTol',2.5e-13,'AbsTol',1e-22); %set tolerences
%[t,y]=ode113(@f1,tspan,y0,options)
y=yg;
t=tg;
plot(y(:,1),y(:,2),'k');
title('Scaled Ganymede Rotating Frame');
[m,p]=size(y);
fg=y(m,:);
p=y;


rg=1.070e6;
vg=10.909;
wg=1.02e-5;
Tg=7.15*24*60*60;
re=6.711e5;
ve=13.780;
we=2.047e-5;
Te=3.55*24*60*60;

%subscript "r" denotes unscaled rotating frame values
pr(:,1)=y(:,1)*rg;
pr(:,2)=y(:,2)*rg;
pr(:,3)=0*rg;
vr(:,1)=y(:,3)*vg;
vr(:,2)=y(:,4)*vg;
vr(:,3)=0*vg;
tr=(Tg/2*pi)*t;

figure
%ganymede to inertial
plot(0,0,'k*');
hold on;
plot(pr(:,1),pr(:,2),'r');
title(' Unscaled Ganymede Rotating Frame');


t=tr;
%Ganymede rotating frame to inertial
for i=1:numSteps
theetag(i,1)=(wg*t(i));
TG2I=[cos(theetag(i)) -sin(theetag(i)) 0;sin(theetag(i)) cos(theetag(i)) 0;0 0 1];
Pgi(i,:)=TG2I*(pr(i,:))';
end

o1=Pgi;

figure
plot(Pgi(:,1),Pgi(:,2),'b');
title('Inertial Frame');
hold on;


%scaled inertial frame
pi1(:,1)=Pgi(:,1)/rg;
pi1(:,2)=Pgi(:,2)/rg;
pi1(:,3)=0*rg;
%vi1(:,1)=Pgi(:,3)/vg;
%vi1(:,2)=Pgi(:,4)/vg;
%vi1(:,3)=0*vg;
ti1=t/(Tg/2*pi);

figure
plot(pi1(:,1),pi1(:,2),'b');
title(' Scaled Inertial Frame');
hold on;


%inertial ganymede to Ganymede rotating frame
for i=1:numSteps
theetag(i,1)=(wg*t(i));
c=Pgi(i,:)';
TI2G=[cos(theetag(i)) sin(theetag(i)) 0;-sin(theetag(i)) cos(theetag(i)) 0;0 0 1];
Pgg(i,:)=TI2G*c;
end
o3=Pgg;




%Inertial to Europa rotating frame

for i=1:numSteps
theetag(i,1)=(we*t(i));
c=Pgi(i,:)';
TI2E=[cos(theetag(i)) sin(theetag(i)) 0;-sin(theetag(i)) cos(theetag(i)) 0;0 0 1];
Pge(i,:)=TI2E*c;
end
o4=Pge;

%scaled europa frame

pr1(:,1)=Pge(:,1)/re;
pr1(:,2)=Pge(:,2)/re;
pr1(:,3)=0*re;
%vr1(:,1)=Pge(:,3)/ve;
%vr1(:,2)=Pge(:,4)/ve;
%vr1(:,3)=0*ve;
tr1=t/(Tg/2*pi);


figure
plot(0,0,'k*');
hold on;
plot(Pgg(:,1),Pgg(:,2),'r');
title('Ganymede Rotating Frame');

figure
plot(0,0,'k*');
hold on;
plot(Pge(:,1),Pge(:,2),'b');
title('Europa Rotating Frame');

figure
plot(0,0,'k*');
hold on;
plot(pr1(:,1),pr1(:,2),'b');
title('Europa Scaled Rotating Frame');


yi=Pgi;
ye=pr1;

end