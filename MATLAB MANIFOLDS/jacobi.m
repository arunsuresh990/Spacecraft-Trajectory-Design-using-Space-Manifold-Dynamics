function c= jacobi(x0,y0,xdot0,ydot0,mu)
%x=input('Enter value of x cordinate:');
%y=input('Enter value of y cordinate:');
%xdot=input('Enter value of xdot :');
%ydot=input('Enter value of ydot :');
%C=input('Enter value of Jacobi Constant :');
x=x0;
y=y0;
xdot=xdot0;
ydot=ydot0;
r1=sqrt((mu+x)^2+(y)^2);
r2=sqrt((mu+x-1)^2+(y)^2);
%Compute the Jacobi Energy
C=x^2+y^2+2*(1-mu)/r1 + 2*mu/r2-(xdot^2+ydot^2)
c=C;
%ydot=sqrt(x^2+y^2+2*(1-mu)/r1 + 2*mu/r2-(xdot^2+C))