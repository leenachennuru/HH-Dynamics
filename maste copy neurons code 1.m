% define dv/dt,dm/dt,dh/dt,dn/dt
KV = @(I,V,m,h,n,GNa,GK,GL,VNa,VK,VL,c) (I - (GNa*(m^3)*h*(V-VNa) + GK*(n^4)*(V-VK)+GL*(V-VL)))/c;
km = @(V,m)  ((0.1*(25-V)/(exp((25-V)/10) - 1))*(1 - m) - m*(4*exp(-V/18)));
kh = @(V,h)  (0.07*exp(-V/20))*(1-h) - (1/(exp((30-V)/10) + 1))*h;
kn = @(V,n)  (kn = ((0.01*(10-V)/(exp((10-V)/10) - 1))*(1-n) -n*(0.125*exp(-V/80)));
t = 50;
dt = 0.01;
imax = t/dt;
%Initialize all matrices to zero
n = input('Enter the number of neurons in the system')
V = cell(1,n);  %Volatage
m = cell(1,n);  %m var
h = cell(1,n);  %h var
n = cell(1,n);  %n var
I = cell(1,n);  %External current
for i = 1:n
    V{1,i} = zeros(1,imax);
    m{1,i} = zeros(1,imax);
    h{1,i} = zeros(1,imax);
    n{1,i} = zeros(1,imax);
    I(1,i) = zeros(1,imax);
end
%Fix the value of parameters of H-H equations & values of simulation variables & Initial state values of  V
GNa=120 ;
GK = 36;
GL=0.3;
R = 10;
VNa = 115;
VK = -12;
VL=10.5995;
c = 1;
I{1,1}(1) = 0; % External current into neuron 1
V{1,1}(1) = 0; % Initial Voltage value at t = 0

icmax = 1;     % # different current (external) values for which trajectories are computed
Istep = 20;    % Step size between successive current (external) values
ivmax = 1;     % # different initial V's that we want to study
Vstep = 0;     % Step size between succeesive V-values

% Fix the initial values of m, h, n for neuron 1 & neuron 2 (taking initial dm/dt=0, dh/dt=0, dn/dt=0
alpham = [];
betam = [];
alphah = [];
betah = [];
alphan = [];
betan = [];
for i = 1:n
    alpham(i) = (0.1*(25-V{1,i}(1)))/(exp((25-V{1,i}(1))/10) - 1);
    betam(i) =  4*exp(-V{1,i}(1)/18);
    aplhah(i) = 0.07*exp(-V{1,i}(1)/20);
    betah(i) = 1/(exp((30-V{1,i}(1))/10) + 1);
    