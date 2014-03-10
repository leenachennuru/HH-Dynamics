clear all
% define dv/dt,dm/dt,dh/dt,dn/dt
KV = @(I,V,m,h,n,GNa,GK,GL,VNa,VK,VL,c) ((I - ((GNa*(m^3)*h*(V-VNa)) + (GK*(n^4)*(V-VK))+(GL*(V-VL))))/c);
km = @(V,m)  ((0.1*(25-V)/(exp((25-V)/10) - 1))*(1 - m) - m*(4*exp(-V/18)));
kh = @(V,h)  (0.07*exp(-V/20))*(1-h) - (1/(exp((30-V)/10) + 1))*h;
kn = @(V,n)  ((0.01*(10-V)/(exp((10-V)/10) - 1))*(1-n) -n*(0.125*exp(-V/80)));
p = 10;
dt = 0.001;
imax = p/dt;

%Initialize all matrices to zero
t = input('Enter the number of neurons in the system');
VV = cell(1,t);  %Volatage
mm = cell(1,t);  %m var
alpham = cell(1,t);
betam = cell(1,t);
manal = cell(1,t);
m2anal = cell(1,t);

h = 0.596; 
n = 0.318;
II = cell(1,t);  %External current

krg = cell(4,t); %Runge Kutta computations for all neurons
mrg = cell(4,t);

Vsave = cell(1,t);
msave = cell(1,t);

Jion = cell(1,t);
JNa  = cell(1,t);
JK   = cell(1,t);
JL   = cell(1,t);
Netcurrent = cell(1,t);

for i = 1:t
    VV{1,i}     = zeros(1,imax);
    mm{1,i}     = zeros(1,imax);
    alpham{1,i} = zeros(1,imax);
    betam{1,i}  = zeros(1,imax);
   
    manal{1,i}   = zeros(1,imax);
    m2anal{1,i}   = zeros(1,imax);
end

%Fix the value of parameters of H-H equations & values of simulation variables & Initial state values of  V
GNa = 120 ;
GK = 36;
GL=0.3;
R = 10;
VNa = 115;
VK = -12;
VL=10.5995;
c = 1;

VV{1,1}(1) = 0; % Initial Voltage value at t = 0

icmax = 1;     % # different current (external) values for which trajectories are computed
Istep = 20;    % Step size between successive current (external) values
ivmax = 1;     % # different initial V's that we want to study
Vstep = 0;     % Step size between succeesive V-values

% Fix the initial values of m, h, n for neuron 1 & neuron 2 (taking initial dm/dt=0, dh/dt=0, dn/dt=0

for i = 1:imax  %Time step             
    for k = 1:t      %Updating all the neurons
        
        Iext = 20;
        
        krg{1,k}(i) = dt*(KV((Iext),VV{1,k}(i),mm{1,k}(i),h,n,GNa,GK,GL,VNa,VK,VL,c))  ;      % k refers to the neurons   
        mrg{1,k}(i) = dt*(km(VV{1,k}(i),mm{1,k}(i)));     % i refers to the number of time steps
        
        
        krg{2,k}(i) = dt*(KV(Iext ,VV{1,k}(i)+(0.5*krg{1,k}(i)),mm{1,k}(i)+(0.5*mrg{1,k}(i)),h,n,GNa,GK,GL,VNa,VK,VL,c));
        mrg{2,k}(i) = dt*(km(VV{1,k}(i)+(0.5*krg{1,k}(i)),mm{1,k}(i)+(0.5*mrg{1,k}(i))));
       
        
        krg{3,k}(i) = dt*(KV(Iext ,VV{1,k}(i)+(0.5*krg{2,k}(i)),mm{1,k}(i)+(0.5*mrg{2,k}(i)),h,n,GNa,GK,GL,VNa,VK,VL,c));
        mrg{3,k}(i) = dt*(km((VV{1,k}(i)+(0.5*krg{2,k}(i))),(mm{1,k}(i)+(0.5*mrg{2,k}(i)))));
        
        
        krg{4,k}(i) = dt*(KV(Iext ,VV{1,k}(i)+krg{3,k}(i),mm{1,k}(i)+mrg{3,k}(i),h,n,GNa,GK,GL,VNa,VK,VL,c));
        mrg{4,k}(i) = dt*(km(VV{1,k}(i)+krg{3,k}(i),mm{1,k}(i)+mrg{3,k}(i)));
        
        
        VV{1,k}(i+1) = VV{1,k}(i) + ((krg{1,k}(i) + 2*krg{2,k}(i) + 2*krg{3,k}(i) + krg{4,k}(i))/6);
        mm{1,k}(i+1) = mm{1,k}(i) + ((mrg{1,k}(i) + 2*mrg{2,k}(i) + 2*mrg{3,k}(i) + mrg{4,k}(i))/6);
        
        
        alpham{1,k}(i) = (0.1*(25-VV{1,k}(i)))/(exp((25-VV{1,k}(i))/10) - 1);
        betam{1,k}(i) =  4*exp(-VV{1,k}(i)/18);
    end
    
    
    for b = 1:t
        manal{1,b}(i+1) = ((Iext-(GK*(n^4))*(VV{1,b}(i+1)-VK))-(GL*(VV{1,b}(i+1))))/(GNa*h*(VV{1,b}(i+1)-VNa));
        m2anal{1,b}(i+1)= (alpham{1,b}(i)/(alpham{1,b}(i) + betam{1,b}(i)));
    end
end
figure,plot(VV{1,1},manal{1,1},VV{1,1},m2anal{1,1})
        