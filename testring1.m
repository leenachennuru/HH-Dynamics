clear all
clc
% define dv/dt,dm/dt,dh/dt,dn/dt
KV = @(I,V,m,h,n,GNa,GK,GL,VNa,VK,VL,c) ((I - ((GNa*(m^3)*h*(V-VNa)) + (GK*(n^4)*(V-VK))+(GL*(V-VL))))/c);
km = @(V,m)  ((((0.1*(25-V))/(exp((25-V)/10) - 1))*(1 - m)) - (m*(4*exp(-V/18))));
kh = @(V,h)  ((0.07*exp(-V/20))*(1-h)) - ((1/(exp((30-V)/10) + 1))*h);
kn = @(V,n)  ((0.01*(10-V)/(exp((10-V)/10) - 1))*(1-n) -n*(0.125*exp(-V/80)));
p = 100;
dt = 0.01;
imax = p/dt;
Iext = 50;     %External current

%Initialize all matrices to zero
ring = input('Do you want to connect the neurons in a ring? If yes enter 0 . If no enter 1')

if (ring == 0)
t = input('Enter the number of neurons in the system')
VV = cell(1,t);  %Volatage
mm = cell(1,t);  %m var
hh = cell(1,t);
nn = cell(1,t);
alpham = cell(1,t);
betam  = cell(1,t);
alphah = cell(1,t);
betah  = cell(1,t);
alphan = cell(1,t);
betan  = cell(1,t);
manal  = cell(1,t);
m2anal = cell(1,t);
II     = zeros(1,t);
       
  
krg = cell(4,t); %Runge Kutta computations for all neurons
mrg = cell(4,t);
hrg = cell(4,t);
nrg = cell(4,t);

Vsave = cell(1,t);
msave = cell(1,t);
hsave = cell(1,t);
nsave = cell(1,t);

Jion = cell(1,t);
JNa  = cell(1,t);
JK   = cell(1,t);
JL   = cell(1,t);
Netcurrent = cell(1,t);

for i = 1:t
    VV{1,i}     = zeros(1,imax);
    mm{1,i}     = zeros(1,imax);
    mm{1,i}(1)  = 0.053;
    hh{1,i}     = zeros(1,imax);
    hh{1,i}(1)  = 0.596;
    nn{1,i}     = zeros(1,imax);
    nn{1,i}(1)  = 0.318;
    
    alpham{1,i} = zeros(1,imax);
    betam{1,i}  = zeros(1,imax);
    alphah{1,i} = zeros(1,imax);
    betah{1,i}  = zeros(1,imax);
    alphan{1,i} = zeros(1,imax);
    betan{1,i}  = zeros(1,imax);
   
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
        
        %for s = 4:1:t-2
        %II(s) = (-2*VV{1,s}(i)+(VV{1,s-2}(i)+ VV{1,s+2}(i)))/R;  % uncomment this in order to add as many neurons as u need and comment line number 100
        %end
        
        II(2) = (-5*VV{1,2}(i) + (VV{1,1}(i) + VV{1,4}(i)+ VV{1,3}(i)+VV{1,5}(i)+VV{1,6}(i)))/R;
        II(3) = (-5*VV{1,3}(i) + (VV{1,1}(i) + VV{1,5}(i)+ VV{1,2}(i)+VV{1,4}(i)+VV{1,6}(i)))/R;
        II(1) = Iext - (((5*VV{1,1}(i)) - (VV{1,2}(i)+ VV{1,3}(i)+ VV{1,4}(i)+VV{1,5}(i)+VV{1,6}(i)))/R);
        II(5) = (-5*VV{1,5}(i) + (VV{1,6}(i)+ VV{1,3}(i)+ VV{1,2}(i) + VV{1,1}(i)+VV{1,4}(i)))/R;
        II(6) = (-5*VV{1,6}(i) + (VV{1,5}(i) + VV{1,4}(i)+ VV{1,3}(i)+ VV{1,2}(i)+VV{1,1}(i)))/R;
        II(4) = (-5*VV{1,4}(i) + (VV{1,2}(i) + VV{1,6}(i)+ VV{1,3}(i)+ VV{1,1}(i)+ VV{1,5}(i)))/R;
         
        krg{1,k}(i) = dt*(KV((II(k)),VV{1,k}(i),mm{1,k}(i),hh{1,k}(i),nn{1,k}(i),GNa,GK,GL,VNa,VK,VL,c))  ;      % k refers to the neurons   
        mrg{1,k}(i) = dt*(km(VV{1,k}(i),mm{1,k}(i)));     % i refers to the number of time steps
        hrg{1,k}(i) = dt*(kh(VV{1,k}(i),hh{1,k}(i))); 
        nrg{1,k}(i) = dt*(kn(VV{1,k}(i),nn{1,k}(i))); 
        
        
        krg{2,k}(i) = dt*(KV(II(k) ,VV{1,k}(i)+(0.5*krg{1,k}(i)),mm{1,k}(i)+(0.5*mrg{1,k}(i)),hh{1,k}(i)+(0.5*hrg{1,k}(i)),nn{1,k}(i)+(0.5*nrg{1,k}(i)),GNa,GK,GL,VNa,VK,VL,c));
        mrg{2,k}(i) = dt*(km(VV{1,k}(i)+(0.5*krg{1,k}(i)),mm{1,k}(i)+(0.5*mrg{1,k}(i))));
        hrg{2,k}(i) = dt*(kh(VV{1,k}(i)+(0.5*krg{1,k}(i)),hh{1,k}(i)+(0.5*hrg{1,k}(i))));
        nrg{2,k}(i) = dt*(kn(VV{1,k}(i)+(0.5*krg{1,k}(i)),nn{1,k}(i)+(0.5*nrg{1,k}(i))));
       
        
        krg{3,k}(i) = dt*(KV(II(k) ,VV{1,k}(i)+(0.5*krg{2,k}(i)),mm{1,k}(i)+(0.5*mrg{2,k}(i)),hh{1,k}(i)+(0.5*hrg{2,k}(i)),nn{1,k}(i)+(0.5*nrg{2,k}(i)),GNa,GK,GL,VNa,VK,VL,c));
        mrg{3,k}(i) = dt*(km((VV{1,k}(i)+(0.5*krg{2,k}(i))),(mm{1,k}(i)+(0.5*mrg{2,k}(i)))));
        hrg{3,k}(i) = dt*(kh((VV{1,k}(i)+(0.5*krg{2,k}(i))),(hh{1,k}(i)+(0.5*hrg{2,k}(i)))));
        nrg{3,k}(i) = dt*(kn((VV{1,k}(i)+(0.5*krg{2,k}(i))),(nn{1,k}(i)+(0.5*nrg{2,k}(i)))));
        
        
        krg{4,k}(i) = dt*(KV(II(k) ,VV{1,k}(i)+krg{3,k}(i),mm{1,k}(i)+mrg{3,k}(i),hh{1,k}(i)+hrg{3,k}(i),nn{1,k}(i)+nrg{3,k}(i),GNa,GK,GL,VNa,VK,VL,c));
        mrg{4,k}(i) = dt*(km(VV{1,k}(i)+krg{3,k}(i),mm{1,k}(i)+mrg{3,k}(i)));
        hrg{4,k}(i) = dt*(kh(VV{1,k}(i)+krg{3,k}(i),hh{1,k}(i)+hrg{3,k}(i)));
        nrg{4,k}(i) = dt*(kn(VV{1,k}(i)+krg{3,k}(i),nn{1,k}(i)+nrg{3,k}(i)));
        
        
        VV{1,k}(i+1) = VV{1,k}(i) + ((krg{1,k}(i) + 2*krg{2,k}(i) + 2*krg{3,k}(i) + krg{4,k}(i))/6);
        mm{1,k}(i+1) = mm{1,k}(i) + ((mrg{1,k}(i) + 2*mrg{2,k}(i) + 2*mrg{3,k}(i) + mrg{4,k}(i))/6);
        hh{1,k}(i+1) = hh{1,k}(i) + ((hrg{1,k}(i) + 2*hrg{2,k}(i) + 2*hrg{3,k}(i) + hrg{4,k}(i))/6);
        nn{1,k}(i+1) = nn{1,k}(i) + ((nrg{1,k}(i) + 2*nrg{2,k}(i) + 2*nrg{3,k}(i) + nrg{4,k}(i))/6);
        
        
        alpham{1,k}(i) = (0.1*(25-VV{1,k}(i)))/(exp((25-VV{1,k}(i))/10) - 1);
        betam{1,k}(i) =  4*exp(-VV{1,k}(i)/18);
        alphah{1,k}(i) = 0.07*exp(-VV{1,k}(i)/20);
        betah{1,k}(i) =  1/(exp((30-VV{1,k}(i))/10) + 1);
        alphan{1,k}(i) = 0.01*(10-VV{1,k}(i))/(exp((10-VV{1,k}(i))/10) - 1);
        betan{1,k}(i) =  0.125*exp(-VV{1,k}(i)/80);
        
        Jion{1,k}(i) = GNa*((mm{1,k}(i))^3)*hh{1,k}(i)*(VV{1,k}(i)-VNa) + GK*((nn{1,k}(i))^4)*(VV{1,k}(i)-VK);
        JNa{1,k}(i) = GNa*((mm{1,k}(i))^3)*hh{1,k}(i)*(VV{1,k}(i)-VNa);
        JK{1,k}(i) = GK*((nn{1,k}(i))^4)*(VV{1,k}(i)-VK);
        JL{1,k}(i) = GL*(VV{1,k}(i)-VL);
        
    end 
end

else if (ring == 1)
       fprintf('Please enter the number of neurons')
       
    else
        fprintf('you have entered an invalid input. Please choose either 0 or 1')
    end
end

figure,subplot(2,3,1); plot(VV{1,1})
subplot(2,3,2);plot(VV{1,2})
subplot(2,3,3);plot(VV{1,3})
subplot(2,3,4);plot(VV{1,4})
subplot(2,3,5);plot(VV{1,5})
subplot(2,3,6);plot(VV{1,6})

%figure,subplot(2,3,1);plot(VV{1,7})
%subplot(2,3,2);plot(VV{1,8})
%subplot(2,3,3);plot(VV{1,9})
%subplot(2,3,4);plot(VV{1,10})
%subplot(2,3,5);plot(VV{1,11})
%subplot(2,3,6);plot(VV{1,12})