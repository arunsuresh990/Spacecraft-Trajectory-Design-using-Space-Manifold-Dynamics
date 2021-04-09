function hetroclinicsplittrajectory1
format long;
mu = 7.802e-5;
mu1=1-mu;
epsilon=2*10^(-4);

%i0=[ 0.99992198;.01557542;- 0.056623322622419;-.00037135];
y0=[1.034731870033446 ;  0.028798141703113;  -0.009107415607789 ; -0.013919377706517]
c= jacobi(y0(1),y0(2),y0(3),y0(4),mu)
hillscurve(c);
numSteps=200;
p=1
ti=0;
tend=6;
k=tend/numSteps

for i=0:k:tend-k
x(p,:)=y0;   
tf=ti+k;
tspan = [ti tf];
options = odeset('RelTol',1e-12,'AbsTol',1e-8);             
[t,y] = ode45(@f,tspan,y0,options);
[m,z]=size(y);
y0=y(m,:);
p=p+1;
plot (mu1,0,'*');
ti=tf;
end
plot(x(:,1),x(:,2),'m');
function dydt = f(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(2))^2);
r2=sqrt((-1+mu+y(1))^2+(y(2))^2);
dydt=[y(3);
      y(4);
      (y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      (y(2)-2*y(3)-(1-mu)*y(2)/(r1)^3-mu*y(2)/(r2)^3)];
end
end