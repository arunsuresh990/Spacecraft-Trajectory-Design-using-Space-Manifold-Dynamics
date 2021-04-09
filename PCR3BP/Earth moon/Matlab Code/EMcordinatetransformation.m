function [yg]= EMcordinatetransformation(ye,te)
%given a vector Xe= [x y xdot ydot] and its time in JE frame, calculates its cordinates in JG frame Xg 

Qe=ye';
pe=te;%present time of ganymede
Tg=2.361e6;%
Te=3.147e7;%
Teg=Te/Tg;
Lg=3.850e5;
Le=1.496e8;
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
dg=[1-mug;0;0;0];
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