clc
clear all
kx = @(x,y) (y+2);
ky = @(x,y) (((-(exp((-((sqrt(x^2+y^2))-1)^2)/10)+ exp(-(x^2+ y^2)/(10))+((exp((-((sqrt(x^2 + y^2))-2)^2)/10))*(-y)))*(x^2 -1/3)*y))- x^5);
p = 10;
dt = 0.01;
imax = p/dt;


%Initialize all matrices to zero

t = input('Enter the number of neurons in the system');
xx = cell(1,t);  %Volatage
yy = cell(1,t);  %m var

       
  
xrg = cell(4,t); %Runge Kutta computations for all neurons
yrg = cell(4,t);


xsave = cell(1,t);
ysave = cell(1,t);



for i = 1:t
    xx{1,i}     = zeros(1,imax);
    yy{1,i}     = zeros(1,imax);   
end

xx{1,1}(1) = 10;
yy{1,1}(1) = 10;

for i = 1:imax  %Time step             
    for k = 1:t      %Updating all the neurons

                 
        xrg{1,k}(i) = dt*(kx(xx{1,k}(i),yy{1,k}(i)));     % i refers to the number of time steps
        yrg{1,k}(i) = dt*(ky(xx{1,k}(i),yy{1,k}(i)));     % k refers to the neurons
        
        
        
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



figure, plot(xx{1,1})
figure,plot(yy{1,1})
%subplot(2,2,3);plot(VV{1,3})
%subplot(2,2,4);plot(VV{1,4})
%subplot(2,3,5);plot(VV{1,5})
%subplot(2,3,6);plot(VV{1,6})

%figure,subplot(2,3,1);plot(VV{1,7})
%subplot(2,3,2);plot(VV{1,8})
%subplot(2,3,3);plot(VV{1,9})
%subplot(2,3,4);plot(VV{1,10})
%subplot(2,3,5);plot(VV{1,11})
%subplot(2,3,6);plot(VV{1,12})