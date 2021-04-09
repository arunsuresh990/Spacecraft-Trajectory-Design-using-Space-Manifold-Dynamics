function [ye]= JEcordinatetransformation(yg,te)
%given a vector Xg= [x y xdot ydot] and its time in JG frame, calculates its cordinates in JE frame Xe 

Qg=yg';
pg=te;%present time in ganymede frame
Tg=6.165e5;%(7.15*24*60*60)ganymede unscaled time period
Te=3.060e5;%(3.55*24*60*60)europa unscaled time period
Tge=Tg/Te;
Lg=1.07e6;
Le=6.711e5;
Lge=Lg/Le;


%l=pg*Tg/(2*pi);%ganymede unscaled time
%omegaG=1.02e-5;
theetaG=pg;
c=cos(theetaG);
s=sin(theetaG);
A=[c -s;s c];
B=[0 0;0 0];
C=[-s -c;c -s];
D=[c -s;s c];
R=[A B;C D];
P=inv(R);
mug= 7.802e-5;
dg=[-mug;0;0;0];
mue =2.523e-5;
de=[-mue;0;0;0];
kg=Qg-dg;
Qgi=R*kg;


Qei(1)=Lge*Qgi(1);
Qei(2)=Lge*Qgi(2);
Qei(3)=(Lge/Tge)*Qgi(3);
Qei(4)=(Lge/Tge)*Qgi(4);
u=Qei';
pe=Tge*pg; 
%j=pe*Te/(2*pi);%europa unscaled time


ke=Qei';
Qe=(P*ke)+de;
%ye=Qgi';
ye=Qe';