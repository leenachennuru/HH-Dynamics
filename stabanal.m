clc
clear all

V = 0:0.01:115;
m = 0:1/11500:1;
h = 0:1/11500:1;
n = 0:1/11500:1;

VVm = [V;m];
VVh = [V;h];
VVn = [V;n];

GNa = 120 ;
GK = 36;
GL=0.3;
R = 10;
VNa = 115;
VK = -12;
VL=10.5995;
c = 1;
Iext = 400;

am = @(V)  ((0.1*(25-V))/(exp((25-V)/10) - 1));
bm = @(V)  (4*exp(-V/18)) ;
ah = @(V)  (0.07*exp(-V/20));
bh = @(V)  (1/(exp((30-V)/10) + 1));
an = @(V)  (0.01*(10-V)/(exp((10-V)/10) - 1));
bn = @(V)  (0.125*exp(-V/80));
mm = @(V)  ((0.1*(25-V))/(exp((25-V)/10) - 1))/ (((0.1*(25-V))/(exp((25-V)/10) - 1)) +  ( 4*exp(-V/18)));
hh = @(V)  (0.07*exp(-V/20))/((0.07*exp(-V/20))+ ( 1/(exp((30-V)/10) + 1)));
nn = @(V)  (0.01*(10-V)/(exp((10-V)/10) - 1))/ ((0.01*(10-V)/(exp((10-V)/10) - 1)) + (0.125*exp(-V/80)));
%mm = @(V)  (am(V))/ ((am(V))+ bm(V));
KVm = @(V,m) ((Iext - ((GNa*(m^3)*hh(V)*(V-VNa)) + (GK*(nn(V)^4)*(V-VK))+(GL*(V-VL))))/c);
KVh = @(V,h) ((Iext - ((GNa*(mm(V)^3)*h*(V-VNa)) + (GK*(nn(V)^4)*(V-VK))+(GL*(V-VL))))/c);
KVn = @(V,n) ((Iext - ((GNa*(mm(V)^3)*hh(V)*(V-VNa)) + (GK*(n^4)*(V-VK))+(GL*(V-VL))))/c);
km = @(V,m)  ((0.1*(25-V)/(exp((25-V)/10) - 1))*(1 - m) - m*(4*exp(-V/18)));
kh = @(V,h)  (0.07*exp(-V/20))*(1-h) - (1/(exp((30-V)/10) + 1))*h;
kn = @(V,n)  ((0.01*(10-V)/(exp((10-V)/10) - 1))*(1-n) -n*(0.125*exp(-V/80)));
Fm = @(VVm)    ((Iext - ((GNa*(VVm(2)^3)*hh(VVm(1))*(VVm(1)-VNa)) + (GK*(nn(VVm(1))^4)*(VVm(1)-VK))+(GL*(VVm(1)-VL))))/c);
Fh = @(VVh)    ((Iext - ((GNa*(VVh(2)^3)*hh(VVh(1))*(VVh(1)-VNa)) + (GK*(nn(VVh(1))^4)*(VVh(1)-VK))+(GL*(VVh(1)-VL))))/c);
Fn = @(VVn)    ((Iext - ((GNa*(VVn(2)^3)*hh(VVn(1))*(VVn(1)-VNa)) + (GK*(nn(VVn(1))^4)*(VVn(1)-VK))+(GL*(VVn(1)-VL))))/c);
InitialGuess = [18.44;0.3304];
InitialGuessh = [1;1];
InitialGuessn = [1;1];
Options = optimset('Display','iter');
XYm = fsolve(Fm, InitialGuess, Options);
ShouldBeZero = Fm(XYm)

V  = 18.46;   
m  = 0.3304; 
h  = 0.1038;  %correspond to current = 100
n  = 0.5995;

Vm = [V, m];
ShouldBeZerom = Fm(Vm)
