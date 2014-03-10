% define dv/dt,dm/dt,dh/dt,dn/dt
KV = @(I,V,m,h,n,GNa,GK,GL,VNa,VK,VL,c) ((I - ((GNa*(m^3)*h*(V-VNa)) + (GK*(n^4)*(V-VK))+(GL*(V-VL))))/c);
km = @(V,m)  ((0.1*(25-V)/(exp((25-V)/10) - 1))*(1 - m) - m*(4*exp(-V/18)));
kh = @(V,h)  (0.07*exp(-V/20))*(1-h) - (1/(exp((30-V)/10) + 1))*h;
kn = @(V,n)  ((0.01*(10-V)/(exp((10-V)/10) - 1))*(1-n) -n*(0.125*exp(-V/80)));
t = 50;
dt = 0.01;
imax = t/dt;

%Initialize all matrices to zero
n = input('Enter the number of neurons in the system');
VV = cell(1,n);  %Volatage
mm = cell(1,n);  %m var
hh = cell(1,n);  %h var
nn = cell(1,n);  %n var
II = cell(1,n);  %External current

krg = cell(4,n); %Runge Kutta computations for all neurons
mrg = cell(4,n);
hrg = cell(4,n);
nrg = cell(4,n);
Vsave = cell(1,n);
msave = cell(1,n);
hsave = cell(1,n);
nsave = cell(1,n);
Jion = cell(1,n);
JNa  = cell(1,n);
JK   = cell(1,n);
JL   = cell(1,n);
Netcurrent = cell(1,n);

for i = 1:n
    VV{1,i} = zeros(1,imax);
    mm{1,i} = zeros(1,imax);
    hh{1,i} = zeros(1,imax);
    nn{1,i} = zeros(1,imax);
    II{1,i} = zeros(1,imax);
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
II{1,1}(1) = 0; % External current into neuron 1
VV{1,1}(1) = 0; % Initial Voltage value at t = 0

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
    alpham(i) = (0.1*(25-VV{1,i}(1)))/(exp((25-VV{1,i}(1))/10) - 1);
    betam(i) =  4*exp(-VV{1,i}(1)/18);
    alphah(i) = 0.07*exp(-VV{1,i}(1)/20);
    betah(i) = 1/(exp((30-VV{1,i}(1))/10) + 1);
    alphan(i)= 0.01*(10-VV{1,i}(1))/(exp((10-VV{1,i}(1))/10) - 1);
    betan(i) = 0.125*exp(-VV{1,i}(1)/80);
    mm{1,i}(1)= (alpham(i)/(alpham(i) + betam(i)));
    hh{1,i}(1)= (alphah(i)/(alphah(i) + betah(i)));
    nn{1,i}(1) = (alphan(i)/(alphan(i) + betan(i)));
end
for ic = 1:icmax
    II{1,1}(1) = II{1,1}(1) + Istep;   %Loop over different values of external current
    for i = 1:imax  %Time step
        for j = 2:n
            II{1,j}(1) = ((VV{1,j-1}(1)- VV{1,j}(1))/R);
            %I{1,i}(1) = 0 uncoupled case
        end
    %Updating all the neurons
    for k = 1:n
        II{1,n+1} = zeros(1,imax);
        krg{1,k}(i) = dt*(KV((II{1,k}(i) - II{1,k+1}(i)),VV{1,k}(i),mm{1,k}(i),hh{1,k}(i),nn{1,k}(i),GNa,GK,GL,VNa,VK,VL,c));           % k refers to the neurons   
        mrg{1,k}(i) = dt*(km(VV{1,k}(i),mm{1,k}(i)));     % i refers to the number of time steps
        hrg{1,k}(i) = dt*(kh(VV{1,k}(i),hh{1,k}(i)));
        nrg{1,k}(i) = dt*(kn(VV{1,k}(i),nn{1,k}(i)));
        
        krg{2,k}(i) = dt*(KV(II{1,k}(i) - II{1,k+1}(i),VV{1,k}(i)+(0.5*krg{1,k}(i)),mm{1,k}(i)+(0.5*mrg{1,k}(i)),hh{1,k}(i)+(0.5*hrg{1,k}(i)),nn{1,k}(i)+(0.5*nrg{1,k}(i)),GNa,GK,GL,VNa,VK,VL,c));
        mrg{2,k}(i) = dt*(km(VV{1,k}(i)+(0.5*krg{1,k}(i)),mm{1,k}(i)+(0.5*mrg{1,k}(i))));
        hrg{2,k}(i) = dt*(kh(VV{1,k}(i)+(0.5*krg{1,k}(i)),hh{1,k}(i)+(0.5*hrg{1,k}(i))));
        nrg{2,k}(i) = dt*(kn(VV{1,k}(i)+(0.5*krg{1,k}(i)),nn{1,k}(i)+(0.5*nrg{1,k}(i))));
        
        krg{3,k}(i) = dt*(KV(II{1,k}(i) - II{1,k+1}(i),VV{1,k}(i)+(0.5*krg{2,k}(i)),mm{1,k}(i)+(0.5*mrg{2,k}(i)),hh{1,k}(i)+(0.5*hrg{2,k}(i)),nn{1,k}(i)+(0.5*nrg{2,k}(i)),GNa,GK,GL,VNa,VK,VL,c));
        mrg{3,k}(i) = dt*(km((VV{1,k}(i)+(0.5*krg{2,k}(i))),(mm{1,k}(i)+(0.5*mrg{2,k}(i)))));
        hrg{3,k}(i) = dt*(kh(VV{1,k}(i)+(0.5*krg{2,k}(i)),hh{1,k}(i)+(0.5*hrg{2,k}(i))));
        nrg{3,k}(i) = dt*(kn(VV{1,k}(i)+(0.5*krg{2,k}(i)),nn{1,k}(i)+(0.5*nrg{2,k}(i))));
        
        krg{4,k}(i) = dt*(KV(II{1,k}(i) - II{1,k+1}(i),VV{1,k}(i)+krg{3,k}(i),mm{1,k}(i)+mrg{3,k}(i),hh{1,k}(i)+hrg{3,k}(i),nn{1,k}(i)+nrg{3,k}(i),GNa,GK,GL,VNa,VK,VL,c));
        mrg{4,k}(i) = dt*(km(VV{1,k}(i)+krg{3,k}(i),mm{1,k}(i)+mrg{3,k}(i)));
        hrg{4,k}(i) = dt*(kh(VV{1,k}(i)+krg{3,k}(i),hh{1,k}(i)+hrg{3,k}(i)));
        nrg{4,k}(i) = dt*(kn(VV{1,k}(i)+krg{3,k}(i),nn{1,k}(i)+nrg{3,k}(i)));
        
        VV{1,k}(i+1) = VV{1,k}(i) + ((krg{1,k}(i) + 2*krg{2,k}(i) + 2*krg{3,k}(i) + krg{4,k}(i))/6);
        mm{1,k}(i+1) = mm{1,k}(i) + ((mrg{1,k}(i) + 2*mrg{2,k}(i) + 2*mrg{3,k}(i) + mrg{4,k}(i))/6);
        hh{1,k}(i+1) = hh{1,k}(i) + ((hrg{1,k}(i) + 2*hrg{2,k}(i) + 2*hrg{3,k}(i) + hrg{4,k}(i))/6);
        nn{1,k}(i+1) = nn{1,k}(i) + ((nrg{1,k}(i) + 2*nrg{2,k}(i) + 2*nrg{3,k}(i) +nrg{4,k}(i))/6);
        
        Vsave{1,k}(i+1,ic) = VV{1,k}(i+1);
        msave{1,k}(i+1,ic) = VV{1,k}(i+1);
        hsave{1,k}(i+1,ic) = VV{1,k}(i+1);
        nsave{1,k}(i+1,ic) = VV{1,k}(i+1);        
        
        Jion{1,k}(i) = GNa*(((mm{1,k}(i))^3)*hh{1,k}(i)*(VV{1,k}(i)-VNa)) + GK*((nn{1,k}(i))^4)*(VV{1,k}(i)-VK);
        JNa{1,k}(i)  = GNa*((mm{1,k}(i))^3)*hh{1,k}(i)*(VV{1,k}(i)-VNa);
        JK{1,k}(i) = GK*((nn{1,k}(i))^4)*(VV{1,k}(i)-VK);
        JL{1,k}(i) = GL*(VV{1,k}(i)-VL);
    end
    end
end

        
        