function M= monodromy(haloPeriod,i0)
t=haloPeriod;
y0=i0;
%function monodromy
%t=3.161453771661801;
%y0=[ 0.965251237914046;0;0;0.046053055752129];
mu=7.802e-5;
mu1=1-mu;
init1=[y0(1) y0(2) y0(3) y0(4)];
y0 = [y0(1);y0(2);y0(3);y0(4);
    1;0;0;0;
    0;1;0;0;
    0;0;1;0;
    0;0;0;1];
tspan = [0 t];

options = odeset('RelTol',1e-11,'AbsTol',1e-22);
             
[t,y] = ode45(@f,tspan,y0,options);

[m,n]=size(y);
k=y(m,:);
mon=[k(5) k(6) k(7) k(8);k(9) k(10) k(11) k(12);k(13) k(14) k(15) k(16);k(17) k(18) k(19) k(20)];
o=det(mon);
M=mon;

function dydt = f(t,y)
  % Derivative function -- mu and mu1 shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(3))^2);
r2=sqrt((-mu1+y(1))^2+(y(3))^2);
p=1-mu1*(1/r1^3-3*(y(1)+mu)^2/r1^5)-mu*(1/r2^3-3*(y(1)-mu1)^2/r2^5);
r=1-mu1*(1/r1^3-3*(y(3))^2/r1^5)-mu*(1/r2^3-3*(y(3))^2/r2^5);
q=3*mu1*y(3)*(y(1)+mu)/r1^5+3*mu*y(3)*(y(1)-mu1)/r2^5;

dydt=[y(2);
    y(1)+2*y(4)-((1-mu)*(y(1)+mu)/(r1)^3)-(mu*(y(1)-(1-mu))/r2^3);
    y(4);
    y(3)-2*y(2)-((1-mu)*y(3)/(r1)^3)-(mu*y(3)/(r2)^3);
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

