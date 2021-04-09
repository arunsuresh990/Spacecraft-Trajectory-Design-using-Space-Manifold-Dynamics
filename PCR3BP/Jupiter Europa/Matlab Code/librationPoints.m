%LIBERATION POINTS plot 0.836292590899931, 1.156168165905523, -1.005115511606892, 0.487722529000000   0.866025403784439
function [l1,l2,l3,l4,l5] =librationPoints(mu)
%mu = 0.012277471;
l=1-mu;
%L3
p_L3=[1, 2*(mu-l), l^2-4*mu*l+mu^2, 2*mu*l*(l-mu)+(l+mu), mu^2*l^2+2*(mu^2-l^2), l^3+mu^3];
L3roots=roots(p_L3);
%initialize L3 for loop
L3=0;
for i=1:5
if L3roots(i) < -mu
L3=L3roots(i);
end
end
%L1
p_L1=[1, 2*(mu-l), l^2-4*l*mu+mu^2, 2*mu*l*(l-mu)+mu-l, mu^2*l^2+2*(l^2+mu^2), mu^3-l^3];
L1roots=roots(p_L1);
%initialize L1 for loop
L1=0;
for i=1:5
if (L1roots(i) > -mu) && (L1roots(i) < l)
L1=L1roots(i);
end
end
%L2
p_L2=[1, 2*(mu-l), l^2-4*l*mu+mu^2, 2*mu*l*(l-mu)-(mu+l), mu^2*l^2+2*(l^2-mu^2), -(mu^3+l^3)];
L2roots=roots(p_L2);
%initialize L2 for loop
L2=0;
for i=1:5
if (L2roots(i) > -mu) && (L2roots(i) > l)
L2=L2roots(i);
end
end
%L4
L4=[-mu+0.5;
sqrt(3)/2];
%L5
L5=[-mu+0.5;
    
-sqrt(3)/2];
l1=L1';
l2=L2';
l3=L3';
l4=L4';
l5=L5';
end