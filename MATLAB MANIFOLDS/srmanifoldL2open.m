function srmanifoldL2open
mu = 7.802e-5;
mu1=1-mu;
epsilon=2*10^(-4);
haloPeriod= 3.259219766232995;%2.972922596049174 ;%3.161454133275042;
i0=[  1.020104811858536;0;0; 0.053526719998360];
c= jacobi(i0(1),i0(2),i0(3),i0(4));
hillscurve(c);
%i0=[ 0.965251237914046;0;0; 0.046053055752129];% c=3.0075940671432188,halfway  0.971739872159348   0.000000176717469  -0.000000000000000  -0.007702968963519
T=3.1*haloPeriod;
k=40;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);
%t=haloPeriod;
%tspan = [0 t];
options=odeset('RelTol',1e-11,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f,tspan,i0,options);%FUNCTION CALL
plot(y(:,1),y(:,3),'r');
hold on;
mon=monodromy(haloPeriod,i0);%FUNCTION CALL
[V,E] = eigs(mon)
moneig=diag(E);

%get the parameterization points of the periodic orbit
numSteps=300;
h=round((numSteps-1)/k);
for n=1:k
haloPoint(n,1:4)=y((n-1)*h+1,1:4);
haloTimes(n)=t((n-1)*h+1);
end
I=eye(4);
init=haloPoint;
for i=1:40
plot(init(i,1),init(i,3),'k*');
hold on;
end
deltaU=V(1:4,1)
deltaS=V(1:4,4)



%STABLE RIGHT
for i=1:40
xU_p(1,1:4)=init(i,1:4)'+epsilon*I*deltaS;
%P=[0 0 0 0;0 1 0 0;0 0 0 0;0 0 0 1];
%xU_p(1,1:4)=init(1,1:4)'+epsilon*I*deltaU;
y0=xU_p(1,:);
tspan = [0 5];
options = odeset('RelTol',1e-12,'AbsTol',1e-8);             
[t,y] = ode113(@f1,tspan,y0,options);
plot(y(:,1),y(:,3),'b');
title('Restricted three body problem');
ylabel('y(t)');
xlabel('x(t)');
plot (mu,0,'*');
plot (mu1,0,'*');
hold on;
end




 function dydt = f(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(3))^2);
r2=sqrt((-1+mu+y(1))^2+(y(3))^2);
dydt=[y(2);
     (y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      y(4);
     (y(3)-2*y(2)-(1-mu)*y(3)/(r1)^3-mu*y(3)/(r2)^3)];
  end  
function dydt = f1(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(3))^2);
r2=sqrt((-1+mu+y(1))^2+(y(3))^2);
dydt=[-y(2);
     -(y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      -y(4);
     -(y(3)-2*y(2)-(1-mu)*y(3)/(r1)^3-mu*y(3)/(r2)^3)];
  end  
end