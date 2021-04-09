function plotallhillscurve
%Jacobi constants and respective hills curve
mu=1.215e-2;
[l1,l2,l3,l4,l5] =librationPoints(mu)
c1= jacobi(l1,0,0,0);
c2= jacobi(l2,0,0,0);
c3= jacobi(l3,0,0,0);
c4= jacobi(l4(1),l4(2),0,0);
figure
hillscurve(c1);
figure
hillscurve(c2);
figure
hillscurve(c3);
figure
hillscurve(c4)
