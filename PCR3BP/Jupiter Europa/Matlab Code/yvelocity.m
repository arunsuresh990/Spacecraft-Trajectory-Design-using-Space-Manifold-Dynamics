x=input('Enter value of x cordinate:');
y=input('Enter value of y cordinate:');
ydot=input('Enter value of ydot :');
%xdot=input('Enter value of xdot :');
C=input('Enter value of Jacobi Constant :');
mu=input('Enter value of mu :');
%mu=7.802e-5;
r1=sqrt((mu+x)^2+(y)^2);
r2=sqrt((mu+x-1)^2+(y)^2);
%Compute the Jacobi Energy
xdot=sqrt(x^2+y^2+2*(1-mu)/r1 + 2*mu/r2-(ydot^2+C))
%ydot=sqrt(x^2+y^2+2*(1-mu)/r1 + 2*mu/r2-(xdot^2+C))