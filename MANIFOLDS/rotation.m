x=[.4 .3 .2 .1]
y=[0 .04 .09 .16]
plot(x,y,'r-');
hold on;
[m,n]=size(y)
%plot(x(1),x(2),'r*');
t=(60*pi/180);
for i=1:n
    ix(i)=x(i)*cos(t)-y(i)*sin(t);
    iy(i)=x(i)*sin(t)+y(i)*cos(t);
end
p=ix
q=iy
hold on;
plot(ix,iy,'c-');
%axis([-5 5 -5 5]);
