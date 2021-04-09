%LIBERATION POINTS
function lpplotwithmu
mu=7.802e-5;
[l1,l2,l3,l4,l5] =librationPoints(mu);
a=[l1;l2;l3];
plot(a,0,'*');
hold on;
plot(l4(1),l4(2),'*')
hold on;
plot(l5(1),l5(2),'*')
hold on;
M=(1-mu);
J=-mu;
plot(J,0,'ro');
hold on;
plot(M,0,'ro');
xlabel('x in DU');
ylabel('y in DU');
end

