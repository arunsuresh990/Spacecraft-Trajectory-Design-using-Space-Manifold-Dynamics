function hillscurve(c)
%C=input('Enter value of Energy,limits are "3.007640   3.007536   3.0000780   2.999921   2.999921":');
C=c;
mu=7.802e-5;

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
 [c, h] = contour(x, y, z, C,'k.');
 title('Hills Curve of Jupiter Ganymede System');
hold on;
end