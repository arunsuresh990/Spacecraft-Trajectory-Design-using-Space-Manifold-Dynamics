function rotatehillscurve
mu=7.802e-5;
c=3.0058;
[P,Q]=meshgrid(-1.5:0.005:1.5);
r1=((P + mu).^2+Q.^2).^(1/2);
r2=((P + mu-1).^2+Q.^2).^(1/2);
Z=2*[(1/2)*(P.^2 + Q.^2)+(1-mu)./r1+mu./r2];
%[C,h] = contour(1.596*P,1.596*Q,1.596*Z,1.596*c,'k');
[C,h] = contour(P,Q,Z,c,'k');

rotate(get(h,'children'),[0 0 1],-28.6) 
end