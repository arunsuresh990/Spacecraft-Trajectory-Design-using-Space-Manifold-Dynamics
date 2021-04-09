function [yi,ye]= GframetoEframe(yg,tg,num)
numSteps=num;
y=yg;
t=tg;
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

t=tr;
%Ganymede rotating frame to inertial
for i=1:numSteps
theetag(i,1)=(wg*t(i));
TG2I=[cos(theetag(i)) -sin(theetag(i)) 0;sin(theetag(i)) cos(theetag(i)) 0;0 0 1];
Pgi(i,:)=TG2I*(pr(i,:))';
end

o1=Pgi;

%scaled inertial frame
pi1(:,1)=Pgi(:,1)/rg;
pi1(:,2)=Pgi(:,2)/rg;
pi1(:,3)=0*rg;
%vi1(:,1)=Pgi(:,3)/vg;
%vi1(:,2)=Pgi(:,4)/vg;
%vi1(:,3)=0*vg;
ti1=t/(Tg/2*pi);


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

yi=Pgi;
ye=pr1;

end