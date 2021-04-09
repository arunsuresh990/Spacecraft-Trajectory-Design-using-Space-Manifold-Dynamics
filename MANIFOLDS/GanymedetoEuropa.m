function GanymedetoEuropa
format long;
mu = 7.802e-5;
mu1=1-mu;

i01=[ 0.981838015975331   0.746397578   0.18 -0.224673294703547];
p=6.250521544;
i0=EJcordinatetransformation(i01,p)
c1= jacobi(i0(1),i0(2),i0(3),i0(4),mu);
c= 3.00576;
%scaledhillscurve(c);
numSteps=2000;
tspan=linspace(0, 5,numSteps);
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y]=ode113(@f1,tspan,i0,options);
p=y
o=t
[m,z]=size(y);
fg=y(m,:)
for i=1:m
yg=y(i,:);
tg=t(i);
tg1(i,:)=tg;
yg1(i,:)=yg;
ye=JEcordinatetransformation(yg,tg);
Qe(i,:)=ye;
end

%plot(Qe(:,1),Qe(:,2),'r');
hold on;

[r,n]=size(Qe);
p=t;
l=p*3.55*24*60*60/(2*pi);
omega=1.02e-5;
hold on;
for i=1:r
    theeta=omega*l(i,1);
    gx(i)=Qe(i,1)*cos(theeta)-Qe(i,2)*sin(theeta);
    gy(i)=Qe(i,1)*sin(theeta)+Qe(i,2)*cos(theeta);



hold on;
end
plot(gx,gy,'r');
%ganymedetrajectory;
hold on;
t = linspace(0,2*pi);
plot(cos(t),sin(t),'k*');
hold on
t = linspace(0,2*pi);
plot(1.596*cos(t),1.596*sin(t),'k*');
M=(1-mu);
J=-mu;
plot(J,0,'r*');
hold on;
plot(M,0,'r*');
hold on;
plot(1.596*M,0,'r*');



mu =2.523e-5;
mu1=1-mu;
i0=[ 0.978970256289329   0.669749586599251   0.177315345423028  -0.24256920193660];
c2= jacobi(i0(1),i0(2),i0(3),i0(4),mu);
%i01=[-1.222124794736756  -0.000000000000059  -0.003918596097935   0.355638203395082];
c= 3.0028
%hillscurve(c);
numSteps=2000;
tspan=linspace(0,5,numSteps);
options=odeset('RelTol',1e-13,'AbsTol',1e-22); %set tolerences
[t,y2]=ode113(@f,tspan,i0,options);
%plot(y1(:,1),y1(:,2),'b');
[r,n]=size(y2);
fe=y2(r,:)
p=t;
l=p*3.55*24*60*60/(2*pi);
omega=2.047e-5;
hold on;
for i=1:r
    theeta=omega*l(i,1);
    ex(i)=y2(i,1)*cos(theeta)-y2(i,2)*sin(theeta);
    ey(i)=y2(i,1)*sin(theeta)+y2(i,2)*cos(theeta);



end
plot(ex,ey,'b');
figure
c= 3.0058;
scaledhillscurve(c);
c= 3.0028
hillscurve(c);
%plot(y(:,1),y(:,2),'r');
plot(Qe(:,1),Qe(:,2),'k');
hold on
plot(y2(:,1),y2(:,2),'k');
hold on;
plot(J,0,'r*');
hold on;
plot(M,0,'r*');
hold on;
plot(1.596*M,0,'r*');

function dydt = f(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(2))^2);
r2=sqrt((-1+mu+y(1))^2+(y(2))^2);
dydt=[y(3);
      y(4);
      (y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      (y(2)-2*y(3)-(1-mu)*y(2)/(r1)^3-mu*y(2)/(r2)^3)];
  end  
function dydt = f1(t,y)
  % Derivative function -- mu and mustar shared with the outer function.
r1=sqrt((mu+y(1))^2+(y(2))^2);
r2=sqrt((-1+mu+y(1))^2+(y(2))^2);
dydt=[-y(3);
     - y(4);
      -(y(1)+2*y(4)-(1-mu)*(y(1)+mu)/(r1)^3-mu/r2^3*(y(1)-(1-mu)));
      -(y(2)-2*y(3)-(1-mu)*y(2)/(r1)^3-mu*y(2)/(r2)^3)];
end
end