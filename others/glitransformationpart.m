for i=1:numSteps
theetag(i,1)=(wg*tg(i));
TG2I=[cos(theetag(i)) -sin(theetag(i)) 0;sin(theetag(i)) cos(theetag(i)) 0;0 0 1];
c=pg(i,:)'
pgi(i,:)=TG2I*c;
end
figure
o1=pi
%m=Pig
plot(pgi(:,1),pgi(:,2),'k');
hold on;

B=[0;0;wg];
for i=1:numSteps
    A=[pg(i,1);pg(i,2);0];
    cor(i,:)=1*cross(B,A);
end

vgi=vg+cor;

for i=1:numSteps
theetag(i,1)=(wg*tg(i));
TG2I=[cos(theetag(i)) -sin(theetag(i)) 0;sin(theetag(i)) cos(theetag(i)) 0;0 0 1];
v=vgi(i,:)';
Vi(i,:)=TG2I*v;
end


re=6.711e5;
ve=13.780;
we=2.047e-5;
Te=3.55*24*60*60;

tr=(Te/2*pi)*t;
te=tr;
for i=1:numSteps
theetae(i,1)=(we*te(i));
c=pgi(i,:)';
TI2E=[cos(theetae(i)) sin(theetae(i)) 0;-sin(theetae(i)) cos(theetae(i)) 0;0 0 1];
pe(i,:)=TI2E*c;
end

pee=pe/re
figure
hillscurve(c);
plot(pee(:,1),pee(:,2),'b');
