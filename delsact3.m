clc
clear all
delf1v = @(V,m,n,h) (-120*m^3*h-36*n^4-3/10);
delf2v = @(V,m,n,h) (-1/10/(exp(5/2-1/10*V)-1)*(1-m)+1/10*(5/2-1/10*V)/(exp(5/2-1/10*V)-1)^2*(1-m)*exp(5/2-1/10*V)+2/9*m*exp(-1/18*V));
delf3v = @(V,m,n,h) (-7/2000*exp(-1/20*V)*(1-h)-1/10/(exp(3-1/10*V)+1)^2*h*exp(3-1/10*V));
delf4v = @(V,m,n,h) (-1/100/(exp(1-1/10*V)-1)*(1-n)+1/10*(1/10-1/100*V)/(exp(1-1/10*V)-1)^2*(1-n)*exp(1-1/10*V)+1/640*n*exp(-1/80*V));

delf1m = @(V,m,n,h) (-360*m^2*h*(V-115));
delf2m = @(V,m,n,h) (-(5/2-1/10*V)/(exp(5/2-1/10*V)-1)-4*exp(-1/18*V));
delf3m = @(V,m,n,h) (0);
delf4m = @(V,m,n,h) (0);

delf1h = @(V,m,n,h) (-120*m^3*(V-115));
delf2h = @(V,m,n,h) (0);
delf3h = @(V,m,n,h) (-7/100*exp(-1/20*V)-1/(exp(3-1/10*V)+1));
delf4h = @(V,m,n,h) (0);

delf1n = @(V,m,n,h) (-144*n^3*(V+12));
delf2n = @(V,m,n,h) (0);
delf3n = @(V,m,n,h) (0);
delf4n = @(V,m,n,h) (-(1/10-1/100*V)/(exp(1-1/10*V)-1)-1/8*exp(-1/80*V));
 
V  = 18.46;   
m  = 0.3304; 
h  = 0.1038;  %correspond to current = 100
n  = 0.5995;


%V = 31.3;
%m = 0.6563;       %correspond to Iext = 400 
%h = 0.02664;
%n = 0.741; 

Jacmat = [delf1v(V,m,h,n), delf1m(V,m,h,n), delf1h(V,m,h,n), delf1n(V,m,h,n); delf2v(V,m,h,n), delf2m(V,m,h,n), delf2h(V,m,h,n), delf2n(V,m,h,n); delf3v(V,m,h,n), delf3m(V,m,h,n), delf3h(V,m,h,n), delf3n(V,m,h,n); delf4v(V,m,h,n), delf4m(V,m,h,n), delf4h(V,m,h,n), delf4n(V,m,h,n)];
eigv = eig(Jacmat)