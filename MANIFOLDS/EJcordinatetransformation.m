function [yg]= EJcordinatetransformation(ye,te)
%given a vector Xe= [x y xdot ydot] and its time in JE frame, calculates its cordinates in JG frame Xg 

Qe=ye';
pe=te;%present time of ganymede
Tg=6.165e5;%(7.15*24*60*60)ganymede unscaled time period
Te=3.060e5;%(3.55*24*60*60)europa unscaled time period
Teg=Te/Tg;
Lg=1.07e6;
Le=6.711e5;
Leg=Le/Lg;


%l=pe*Te/(2*pi);%europa unscaled time
%omegaE=1.02e-5;
theetaE=pe;
c=cos(theetaE);
s=sin(theetaE);
A=[c -s;s c];
B=[0 0;0 0];
C=[-s -c;c -s];
D=[c -s;s c];
R=[A B;C D];
P=inv(R);
mue =2.523e-5;
de=[-mue;0;0;0];
mug= 7.802e-5;
dg=[-mug;0;0;0];
ke=Qe-de;
Qei=R*ke;


Qgi(1)=Leg*Qei(1);
Qgi(2)=Leg*Qei(2);
Qgi(3)=(Leg/Teg)*Qei(3);
Qgi(4)=(Leg/Teg)*Qei(4);
u=Qgi';
pg=Teg*pe; 
%j=pg*Tg/(2*pi);%ganymede unscaled time


kg=Qgi';
Qg=(P*kg)+dg;%Qg=(R\kg)+dg;
%yg=Qei';
yg=Qg';