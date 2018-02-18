clear all
home

colorspec = '-r';
t = linspace(0,2*pi,100);
A = -1;
B = 3;
C = 4;
D = 5;
E = 12;

r = A.*exp(cos(t)) - B.*cos(C.*t) + (sin(t./E)).^D;

polar(t,r,colorspec)
hold on