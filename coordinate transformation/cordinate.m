function cordinate
%cordinate transformation of ganymede, ganymede rotating frame to 
%inertial frame and back to ganymede rotating frame and to europa rotating
%frame.In ganymede frame, ganymede is fixed, so it is a point.in inertial 
%frame will give an anticlkwise circle. back to ganymede frame will give 
%back the single point. in europa frame it gives a clockwse rotating circle
%(half times exactly opposite to inertial frame as angular velocity of 
%europa is twice that of ganymede).that is ,in inertial frame let angular velocities
%be Wei and Wgi. so Wei=2*Wgi. in ganymede frame angular velocity of europa Weg=.5*Wei as ganymede is 
%already rotating forward at its angular velocity(half that of ganymede).in europa frame angular velocity of ganymede Wge=-1*Wgi as europa is 
%already rotating forward  at its angular velocity(twice that of ganymede). similarly europa is transformed
rg=1.070e6;
wg=1.02e-5;
re=6.711e5;
we=2.047e-5;
Tg=7.15*24*60*60;
Te=3.55*24*60*60;

%ganymede to inertial
plot(0,0,'k*');
hold on;
Pg=[rg;0];
plot(Pg(1),Pg(2),'r*');
%inertial frame
numSteps=200;
tr=linspace(0,5,numSteps)';
t=(Tg/2*pi)*tr;
for i=1:200
theetag(i,1)=(wg*t(i));
TG2I=[cos(theetag(i)) -sin(theetag(i));sin(theetag(i)) cos(theetag(i))];
Pgi(i,:)=TG2I*Pg;
end
%figure
o1=Pgi
%m=Pig
plot(t(:,1),Pgi(:,1),'k');
figure
plot(Pgi(:,1),Pgi(:,2),'k');
hold on;
numSteps=10;
t=linspace(0,Tg/4,numSteps)';
for i=1:10
theetag(i,1)=(wg*t(i));
TG2I=[cos(theetag(i)) -sin(theetag(i));sin(theetag(i)) cos(theetag(i))];
Pgi1(i,:)=TG2I*Pg;
end
%figure
o=Pgi1;
plot(Pgi1(:,1),Pgi1(:,2),'r*');



%Europa to inertial
plot(0,0,'k*');
hold on;
Pe=[re;0];
plot(Pe(1),Pe(2),'b*');
%inertial frame
numSteps=200;
t=linspace(0,Tg/4,numSteps)';
for i=1:200
theetae(i,1)=(we*t(i));
TE2I=[cos(theetae(i)) -sin(theetae(i));sin(theetae(i)) cos(theetae(i))];
Pei(i,:)=TE2I*Pe;
end
%figure
o2=Pei
plot(Pei(:,1),Pei(:,2),'k');
numSteps=10;
T=3.55*24*60*60;
t=linspace(0,Tg/4,numSteps)';
for i=1:10
theetae(i,1)=(we*t(i));
TE2I=[cos(theetae(i)) -sin(theetae(i));sin(theetae(i)) cos(theetae(i))];
Pei1(i,:)=TE2I*Pe;
end
%figure
o=Pei1;
plot(Pei1(:,1),Pei1(:,2),'b*');
title('Inertial  Frame');



%inertial ganymede to Ganymede rotating frame
numSteps=200;
t=linspace(0,Tg/4,numSteps)';
for i=1:200
theetag(i,1)=(wg*t(i));
c=Pgi(i,:)';
TI2G=[cos(theetag(i)) sin(theetag(i));-sin(theetag(i)) cos(theetag(i))];
Pgg(i,:)=TI2G*c;
end
o3=Pgg

%plot(Pig(:,1),Pig(:,2),'r');


%inertial Ganymede to Europa rotating frame
numSteps=200;
t=linspace(0,Tg/4,numSteps)';
for i=1:200
theetag(i,1)=(we*t(i));
%TG2I=[cos(theetag(i)) -sin(theetag(i));sin(theetag(i)) cos(theetag(i))];
%Pgi(i,:)=TG2I*Pg;
c=Pgi(i,:)';
TI2E=[cos(theetag(i)) sin(theetag(i));-sin(theetag(i)) cos(theetag(i))];
Pge(i,:)=TI2E*c;
end
o4=Pge



%inertial Europa to Europa rotating frame
numSteps=200;
t=linspace(0,Tg/4,numSteps)';
for i=1:200
theetae(i,1)=(we*t(i));
%TG2I=[cos(theetag(i)) -sin(theetag(i));sin(theetag(i)) cos(theetag(i))];
%Pgi(i,:)=TG2I*Pg;
c=Pei(i,:)';
TI2E=[cos(theetae(i)) sin(theetae(i));-sin(theetae(i)) cos(theetae(i))];
Pee(i,:)=TI2E*c;
end
o5=Pee
%figure
%plot(Pig(:,1),Pig(:,2),'r');

%inertial Europa to Ganymede rotating frame
numSteps=200;
t=linspace(0,Tg/4,numSteps)';
for i=1:200
theetae(i,1)=(wg*t(i));
%TG2I=[cos(theetag(i)) -sin(theetag(i));sin(theetag(i)) cos(theetag(i))];
%Pgi(i,:)=TG2I*Pg;
c=Pei(i,:)';
TI2G=[cos(theetae(i)) sin(theetae(i));-sin(theetae(i)) cos(theetae(i))];
Peg(i,:)=TI2G*c;
end
o6=Peg


figure
plot(0,0,'k*');
hold on;
plot(Pgg(:,1),Pgg(:,2),'r*');
hold on
plot(Peg(:,1),Peg(:,2),'b*');
title('Ganymede Rotating Frame');

figure
plot(0,0,'k*');
hold on;
plot(Pge(:,1),Pge(:,2),'r*');
hold on
plot(Pee(:,1),Pee(:,2),'b*');
plot(Pg(1),Pg(2),'r*');
title('Europa Rotating Frame');