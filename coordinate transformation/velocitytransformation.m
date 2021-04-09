function velocitytransformation
clear;
close all;
xf = [5 5 0]';
w = 1;
i = 1;
xrp = [0 0 0]';
dt = 0.01;
for t=0:dt:pi/4
    th = w*t;
    tr = [cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];
    xr = tr*xf;
    vrn = (xr-xrp)/dt;
    omcr = cross([0 0 1]',xr');
    vrc = 0-omcr;
    xrp = xr;
    tt(i) =t;
    xra(i,:) = xr';
    vrna(i,:) = vrn';
    vrca(i,:) = vrc';
    i = i+1;
end

o1=xra
o2=vrna
o3=vrca
o=tt'
figure(1)
plot(xra(:,1),xra(:,2))

figure(2)
plot(tt(2:i-1),vrna(2:i-1,1),'b');
hold on;
plot(tt(2:i-1),vrca(2:i-1,1),'r');

figure(3)
plot(tt(2:i-1),vrna(2:i-1,2),'b');
hold on;
plot(tt(2:i-1),vrca(2:i-1,2),'r');
end