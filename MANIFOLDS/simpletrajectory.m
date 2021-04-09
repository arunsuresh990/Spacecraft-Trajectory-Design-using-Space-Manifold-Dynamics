function simpletrajectory
% Problem parameters
mu = 7.802e-5;
mu1=1-mu;
c= 3.0028;
hillscurve(c);

y0(1)= 0.978332646874192;
y0(2)= 0.009858514255237;
y0(3)=-0.008036026547563 ;
%y0(4)=-0.043013862776183 ;

r1=sqrt((mu+y0(1))^2+(y0(2))^2);
r2=sqrt((-1+mu+y0(1))^2+(y0(2))^2);
%y0(3)=sqrt((y0(1)^2)+(y0(2)^2)+(2*(mu1)/r1)+(2*mu/r2)-(y0(4)^2)-c);
%y0(4)=sqrt((y0(1)^2)+(y0(2)^2)+(2*(mu1)/r1)+(2*mu/r2)-(y0(3)^2)-c);
%y0 = [y0(1);y0(2);y0(3); y0(4)]
y0=[ 0.968385555854192 ;  0.001299026271682 ;  0.001065500396714 ;  0.016691529435333];


tspan = [0 5];

options = odeset('RelTol',1e-12,'AbsTol',1e-8);             
[t,y] = ode45(@f1,tspan,y0,options);
p=y;
plot(y(:,1),y(:,2),'b');
% figure
% %axis([.9 1.05 -.03 .03]); 
% [m,z]=size(y);
% fg=y(m,:)
% for i=1:m
% yg=y(i,:);
% tg=t(i);
% tg1(i,:)=tg;
% yg1(i,:)=yg;
% ye=JEcordinatetransformation(yg,tg);
% Qe(i,:)=ye;
% end
% %plot(Qe(:,1),Qe(:,2),'r');
% hold on;
% 
% 
% p=t;
% l=p*7.15*24*60*60/(2*pi);
% omega=2.047e-5;
% hold on;
% for i=1:m
%     theeta=omega*l(i,1);
%     gx(i,1)=Qe(i,1)*cos(theeta)-Qe(i,2)*sin(theeta);
%     gy(i,1)=Qe(i,1)*sin(theeta)+Qe(i,2)*cos(theeta);
% 
% 
% plot(gx,gy,'b');
% hold on;
% end
% t = linspace(0,2*pi);
% plot(cos(t),sin(t),'k');
% hold on
% t = linspace(0,2*pi);
% plot(1.596*cos(t),1.596*sin(t),'k');
% M=(1-mu);
% J=-mu;
% plot(J,0,'r*');
% hold on;
% plot(M,0,'r*');
% hold on;
% plot(1.596*M,0,'r*');
% 
% title('Restricted three body problem');
% ylabel('y(t)');
% xlabel('x(t)');
% hold on;
% %plot (1.029,0,'*');
% hold on;
% plot (mu1,0,'*');
% 
% plot(Qe(:,1),Qe(:,2),'r');
% c= 3.0058;
% %scaledhillscurve(c);
% plot(1.596*M,0,'r*');  % -----------------------------------------------------------------------
%   % Nested functions -- problem parameters provided by the outer function.
%   %

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
  % -----------------------------------------------------------------------

 

  % -----------------------------------------------------------------------

end  % liaorbitode

  
  