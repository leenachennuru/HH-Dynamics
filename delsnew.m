clear all
clc
GNa = 120 ;
GK = 36;
GL=0.3;
R = 10;
VNa = 115;
VK = -12;
VL= 10.5995;
c = 1;
Iext = 400;

syms V m h n 
f1 = ((Iext - ((GNa*(m^3)*h*(V-VNa)) + (GK*(n^4)*(V-VK))+(GL*(V-VL))))/c);
f2 = ((0.1*(25-V)/(exp((25-V)/10) - 1))*(1 - m) - m*(4*exp(-V/18)));
f3 = (0.07*exp(-V/20))*(1-h) - (1/(exp((30-V)/10) + 1))*h;
f4 = ((0.01*(10-V)/(exp((10-V)/10) - 1))*(1-n) - n*(0.125*exp(-V/80)));

delf1v = diff(f1,V)
delf1m = diff(f1,m)
delf1h = diff(f1,h)
delf1n = diff(f1,n)

delf2v = diff(f2,V)
delf2m = diff(f2,m)
delf2h = diff(f2,h)
delf2n = diff(f2,n)

delf3v = diff(f3,V)
delf3m = diff(f3,m)
delf3h = diff(f3,h)
delf3n = diff(f3,n)

delf4v = diff(f4,V)
delf4m = diff(f4,m)
delf4h = diff(f4,h)
delf4n = diff(f4,n)