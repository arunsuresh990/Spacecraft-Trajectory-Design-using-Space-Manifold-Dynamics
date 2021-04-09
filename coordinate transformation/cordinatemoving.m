function cordinatemoving
%cordinate transformation of ganymede MANIFOLD frame ganymede rotating frame to 
%inertial frame and back to ganymede rotating frame and to europa rotating
%frame.
rg=1.070e6;
wg=1.02e-5;
re=6.711e5;
we=2.047e-5;
Tg=7.15*24*60*60*10;
Te=3.55*24*60*60;



numSteps=200;
t=linspace(0,Tg/4,numSteps)';
r=1;
for i=1:numSteps
Pin(i,1)=2*t(i)+3;%r*cos(t(i));%2*t(i)+3;
Pin(i,2)=3*t(i)+2;%r*sin(t(i));%3*t(i)+2;
%Pin(i,3)=0*tan(t(i));
end
pos=Pin
%yi=Pi(:,2)



%ganymede to inertial
plot(0,0,'k*');
hold on;
%Pg=[rg;0];
%plot(Pg(1),Pg(2),'r*');
plot(Pin(:,1),Pin(:,2),'k');
%inertial frame
numSteps=200;
t=linspace(0,Tg/4,numSteps)';
for i=1:200
theetag(i,1)=(wg*t(i));
TG2I=[cos(theetag(i)) -sin(theetag(i));sin(theetag(i)) cos(theetag(i))];
Pgi(i,:)=TG2I*(Pin(i,:))';
end
%figure
o1=Pgi
%m=Pig
figure
plot(Pgi(:,1),Pgi(:,2),'b');
hold on;
numSteps=10;
t=linspace(0,Tg/4,numSteps)';
for i=1:10
theetag(i,1)=(wg*t(i));
TG2I=[cos(theetag(i)) -sin(theetag(i));sin(theetag(i)) cos(theetag(i))];
Pgi1(i,:)=TG2I*(Pin(i,:))';
end
%figure



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





figure
title('Ganymede Rotating Frame');
plot(0,0,'k*');
hold on;
plot(Pgg(:,1),Pgg(:,2),'r*');


figure
title('Europa Rotating Frame');
plot(0,0,'k*');
hold on;
plot(Pge(:,1),Pge(:,2),'b*');

