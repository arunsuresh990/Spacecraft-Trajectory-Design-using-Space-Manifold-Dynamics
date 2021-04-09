function monodromymatrix%(z,t)
% Problem parameters
mu=1.215e-2;
mu1=1-mu;
delx=1e-10;
hp1=2.972922596049174;
y0(1)= 0.969539190837874 ;
y0(2)=0  ;
y0(3)= 0 ;
y0(4)=0.007530312113471;
c= jacobi(y0(1),y0(2),y0(3),y0(4))

r1=sqrt((mu+y0(1))^2+(y0(2))^2);
r2=sqrt((-mu1+y0(1))^2+(y0(2))^2);
%y0(4)=sqrt((y0(1)^2)+(y0(3)^2)+(2*(mu1)/r1)+(2*mu/r2)-c);
init1=[y0(1) y0(2) y0(3) y0(4)];
y0 = [y0(1);y0(2);y0(3);y0(4);
    1;0;0;0;
    0;1;0;0;
    0;0;1;0;
    0;0;0;1];
tspan = [0 hp1];

options = odeset('RelTol',1e-11,'AbsTol',1e-22);
             
[t,y] = ode45(@f,tspan,y0,options);
plot(y(:,1),y(:,2,'r');
hold on;
plot(y0(1),y0(2),'*');
%plot (.970589422811693,0,'*');
[m,n]=size(y);
k=y(m,:);
fin1=[k(1) k(2) k(3) k(4)];
mon=[k(5) k(6) k(7) k(8);k(9) k(10) k(11) k(12);k(13) k(14) k(15) k(16);k(17) k(18) k(19) k(20)]
gd=det(mon)

[V,D] = eigs(mon);
moneig=diag(D);
%VARIATION
%x(1)=0.965251237914046;
hp2=2.972922596049174;

y0(1)= y0(1)+delx;
y0(2)=y0(2)+delx;
y0(3)=y0(3)+delx;
y0(4)= y0(4)+delx;


%r1=sqrt((mu+y0(1))^2+(y0(3))^2);
%r2=sqrt((-mu1+y0(1))^2+(y0(3))^2);
%y0(4)=sqrt((y0(1)^2)+(y0(3)^2)+(2*(mu1)/r1)+(2*mu/r2)-c);
init2=[y0(1) y0(2) y0(3) y0(4)];
y0 = [y0(1);y0(2);y0(3);y0(4);
    1;0;0;0;
    0;1;0;0;
    0;0;1;0;
    0;0;0;1];
tspan = [0 hp2];

options = odeset('RelTol',1e-11,'AbsTol',1e-22);
             
[t1,y] = ode45(@f,tspan,y0,options);

[m,n]=size(y);
k1=y(m,:);
fin2=[k1(1) k1(2) k1(3) k1(4)];
%mon=[k(5) k(6) k(7) k(7);k(9) k(10) k(11) k(12);k(13) k(14) k(15) k(16);k(17) k(18) k(19) k(20)];
%gd=det(mon);
z=init2-init1;
x=fin2-fin1;
dox0=z.';
doxt=x.'
a=mon*dox0

%x(i+1)=y0(1);
function dydt = f(t,y)
  % Derivative function -- mu and mu1 shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(2))^2);
r2=sqrt((-1+mu+y(1))^2+(y(2))^2);
p=1-mu1*(1/r1^3-3*(y(1)+mu)^2/r1^5)-mu*(1/r2^3-3*(y(1)-mu1)^2/r2^5);
r=1-mu1*(1/r1^3-3*(y(2))^2/r1^5)-mu*(1/r2^3-3*(y(2))^2/r2^5);
q=3*mu1*y(2)*(y(1)+mu)/r1^5+3*mu*y(2)*(y(1)-mu1)/r2^5;

dydt=[y(3);
      y(4);
      (y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      (y(2)-2*y(3)-(1-mu)*y(2)/(r1)^3-mu*y(2)/(r2)^3);
    y(13);
    y(14);
    y(15);
    y(16);
    y(17);
    y(18);
    y(19);
    y(20);
    p*y(5)+q*y(9)+2*y(17);
    p*y(6)+q*y(10)+2*y(18);
    p*y(7)+q*y(11)+2*y(19);
    p*y(8)+q*y(12)+2*y(20);
    q*y(5)+r*y(9)-2*y(13);
    q*y(6)+r*y(10)-2*y(14);
    q*y(7)+r*y(11)-2*y(15);
    q*y(8)+r*y(12)-2*y(16)];
end 
end