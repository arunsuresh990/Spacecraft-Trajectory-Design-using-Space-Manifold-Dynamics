function velocitycordinate
%inertial frame
numSteps=100
t=linspace(0,pi/4,numSteps)'
r=1;
for i=1:numSteps
Pin(i,1)=r*cos(t(i))+4;%r*cos(t(i));%2*t(i)+3;
Pin(i,2)=r*sin(t(i))+5;%r*sin(t(i));%3*t(i)+2;
Pin(i,3)=0*tan(t(i));
end
pos=Pin
%yi=Pi(:,2)

Vxi(1)=0;
Vyi(1)=0;
for i=2:numSteps
Vxi(i)=(Pin(i,1)-Pin(i-1,1))/(t(i)-t(i-1));
Vyi(i)=(Pin(i,2)-Pin(i-1,2))/(t(i)-t(i-1));
Vzi(i)=0;
end
Vi=[Vxi' Vyi' Vzi']



%inertial to rotating frame
wg=1;
theeta(1)=0;
for i=1:numSteps
theeta(i,1)=(wg*t(i));
TI2R=[cos(theeta(i)) sin(theeta(i)) 0;-sin(theeta(i)) cos(theeta(i)) 0;0 0 1];
c=Pin(i,:)';
Pr(i,:)=TI2R*c;
end

for i=1:numSteps
theeta(i,1)=(wg*t(i));
TI2R=[cos(theeta(i)) sin(theeta(i)) 0;-sin(theeta(i)) cos(theeta(i)) 0;0 0 1];
v=Vi(i,:)';
Vir(i,:)=TI2R*v;
end


o3=Pr
o4=Vir


plot(Pin(:,1),Pin(:,2),'b');
hold on
plot(Pr(:,1),Pr(:,2),'r');
axis equal
grid on
Vxr(1)=0;
Vyr(1)=0;
for i=2:numSteps
    Vxr(i)=(Pr(i,1)-Pr(i-1,1))/(t(i)-t(i-1));
    Vyr(i)=(Pr(i,2)-Pr(i-1,2))/(t(i)-t(i-1));
    Vzr(i)=0;
   % V(i)=sqrt(xdot(i)^2+ydot(i)^2)
end

Vr=[Vxr' Vyr' Vzr']
B=[0;0;wg];
for i=1:numSteps
    A=[Pr(i,1);Pr(i,2);0];
    cor(i,:)=-1*cross(B,A);
end
o=cor
cor=Vir+cor


figure
plot(t(2:i-1,1),Vr(2:i-1,1),'r');
hold on
plot(t(2:i-1,1),cor(2:i-1,1),'b');

figure
plot(t(2:i-1,1),Vr(2:i-1,2),'r');
hold on
plot(t(2:i-1,1),cor(2:i-1,2),'b');
