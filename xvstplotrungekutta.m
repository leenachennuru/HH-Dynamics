clc
clear all
% define dv/dt,dm/dt,dh/dt,dn/dt
kx = @(x,y) (y);
ky = @(x,y)  (((-(exp((-((sqrt(x^2+y^2))-1)^2)/10)+exp(-(x^2+ y^2)/(10))+((exp((-((sqrt(x^2 + y^2))-2)^2)/10))*(-y)))*(x^2 -1/3)*y))- x^5);

t = 50;
dt = 0.01;
imax = t/dt;

%Initialize all matrices to zero
n = input('Enter the number of neurons in the system');
xx = cell(1,n);  %Volatage
yy = cell(1,n);  %m var


xrg = cell(4,n); %Runge Kutta computations for all neurons
yrg = cell(4,n);

xsave = cell(1,n);
ysave = cell(1,n);


for i = 1:n
    xx{1,i} = zeros(1,imax);
    yy{1,i} = zeros(1,imax);
    
end


icmax = 1;     % # different current (external) values for which trajectories are computed
Istep = 20;    % Step size between successive current (external) values
ivmax = 1;     % # different initial V's that we want to study
Vstep = 0;     % Step size between succeesive V-values



for i = 1:imax                  %Time step
    for k = 1:n
        xrg{1,k}(i) = dt*(kx(xx{1,k}(i),yy{1,k}(i)));           % k refers to the neurons   
        yrg{1,k}(i) = dt*(ky(xx{1,k}(i),yy{1,k}(i)));     % i refers to the number of time steps
        
        
        xrg{2,k}(i) = dt*(kx(xx{1,k}(i)+(0.5*xrg{1,k}(i)),yy{1,k}(i)+(0.5*yrg{1,k}(i))));
        yrg{2,k}(i) = dt*(ky(xx{1,k}(i)+(0.5*xrg{1,k}(i)),yy{1,k}(i)+(0.5*yrg{1,k}(i))));
        
        
        xrg{3,k}(i) = dt*(kx((xx{1,k}(i)+(0.5*xrg{2,k}(i))),(yy{1,k}(i)+(0.5*yrg{2,k}(i)))));
        yrg{3,k}(i) = dt*(ky((xx{1,k}(i)+(0.5*xrg{2,k}(i))),(yy{1,k}(i)+(0.5*yrg{2,k}(i)))));
        
        
        xrg{4,k}(i) = dt*(kx(xx{1,k}(i)+xrg{3,k}(i),yy{1,k}(i)+yrg{3,k}(i)));
        yrg{4,k}(i) = dt*(ky(xx{1,k}(i)+xrg{3,k}(i),yy{1,k}(i)+yrg{3,k}(i)));
        
        
        xx{1,k}(i+1) = xx{1,k}(i) + ((xrg{1,k}(i) + 2*xrg{2,k}(i) + 2*xrg{3,k}(i) + xrg{4,k}(i))/6);
        yy{1,k}(i+1) = yy{1,k}(i) + ((yrg{1,k}(i) + 2*yrg{2,k}(i) + 2*yrg{3,k}(i) + yrg{4,k}(i))/6);

        
        
      
    end
end
    
figure,plot(xx{1,1})
figure,plot(yy{1,1})

        
        