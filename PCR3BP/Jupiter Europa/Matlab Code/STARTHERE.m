%LIBERATION POINTS
%For a given Mu value, finds all libration points and associaTED JACOBI
%CONSTANTS. then hilkls curves for a user given C is plotted
format long;
mu=input('Enter value of mu:');
[l1,l2,l3,l4,l5] =librationPoints(mu);
L(1)=l1
L(2)=l2
L(3)=l3
L4=l4
L5=l5
for i=1:3
x=L(i);
y=0;
r1=sqrt((mu+x)^2+(y)^2);
r2=sqrt((mu+x-1)^2+(y)^2);
%Compute the Jacobi Energy
C(i)=x^2+y^2+2*(1-mu)/r1 + 2*mu/r2;
end
i=i+1;
x=L4(1);
y=L4(2);
r1=sqrt((mu+x)^2+(y)^2);
r2=sqrt((mu+x-1)^2+(y)^2);
%Compute the Jacobi Energy
C(i)=x^2+y^2+2*(1-mu)/r1 + 2*mu/r2;
C(i+1)=x^2+y^2+2*(1-mu)/r1 + 2*mu/r2

C=input('Enter value of Energy:');

xmin = -1.5;

xmax = +1.5;

ymin = -1.5;

ymax = +1.5;

delta = 0.005;

[x, y] = meshgrid(xmin: delta: xmax, ymin: delta: ymax);

x1 = - mu;

x2 = 1 - mu;

for i = 1: 1: size(x)
    for j = 1: 1: size(x)

        r1sqr = (x(i, j) - x1)^2 + y(i, j)^2;

        r2sqr = (x(i, j) - x2)^2 + y(i, j)^2;

        z(i, j) = (1 - mu) * (r1sqr + 2 / sqrt(r1sqr))+ mu * (r2sqr + 2 / sqrt(r2sqr)) - mu * (1 - mu);
    end
end


lpplotwithmu;
hold on;
 [c, h] = contour(x, y, z, C,'k');
 title('Hills Curve of Jupiter Ganymede System');
 hold on;

