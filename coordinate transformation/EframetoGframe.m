function [yi,yg]= EframetoGframe(ye,te,num)
numSteps=num;
y=ye;
t=te;
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

%subscript "r" denotes unscaled rotating frame values, position velocity
%and time
pr(:,1)=y(:,1)*re;
pr(:,2)=y(:,2)*re;
pr(:,3)=0*re;
vr(:,1)=y(:,3)*ve;
vr(:,2)=y(:,4)*ve;
vr(:,3)=0*ve;
tr=(Te/2*pi)*t;

t=tr;
%Europa rotating frame to inertial ganymede 
for i=1:numSteps
theetae(i,1)=(we*t(i));
TE2I=[cos(theetae(i)) -sin(theetae(i)) 0;sin(theetae(i)) cos(theetae(i)) 0;0 0 1];
Pei(i,:)=TE2I*(pr(i,:))';
end

o1=Pei;

%scaled inertial frame
pi1(:,1)=Pei(:,1)/re;
pi1(:,2)=Pei(:,2)/re;
pi1(:,3)=0*re;
%vi1(:,1)=Pei(:,3)/ve;
%vi1(:,2)=Pei(:,4)/ve;
%vi1(:,3)=0*ve;
ti1=t/(Te/2*pi);


%Inertial to Europa rotating frame
for i=1:numSteps
theetae(i,1)=(we*t(i));
c=Pei(i,:)';
TI2E=[cos(theetae(i)) sin(theetae(i)) 0;-sin(theetae(i)) cos(theetae(i)) 0;0 0 1];
Pee(i,:)=TI2E*c;
end
o3=Pee;

%Inertial to Ganymede rotating frame

for i=1:numSteps
theetae(i,1)=(wg*t(i));
c=Pei(i,:)';
TI2G=[cos(theetae(i)) sin(theetae(i)) 0;-sin(theetae(i)) cos(theetae(i)) 0;0 0 1];
Peg(i,:)=TI2G*c;
end
o4=Peg;


pr1(:,1)=Peg(:,1)/rg;
pr1(:,2)=Peg(:,2)/rg;
pr1(:,3)=0*rg;
%vr1(:,1)=Peg(:,3)/vg;
%vr1(:,2)=Peg(:,4)/vg;
%vr1(:,3)=0*vg;
tr1=t/(Te/2*pi);

yi=Pei;
yg=pr1;

end