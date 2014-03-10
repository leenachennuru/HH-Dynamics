clear all
clc

GNa = 120 ;
GK = 36;
GL=0.3;
R = 10;
VNa = 115;
VK = -12;
VL=10.5995;
c = 1;
Iext = 200;

V1 = 23.8:0.0001:24;
[a,b] = size(V1);
for i = 1:b
    
    alpham1(i) = (0.1*(25-V1(i)))/(exp((25-V1(i))/10) - 1);
    betam1(i) =  4*exp(-V1(i)/18);
    alphah1(i) = 0.07*exp(-V1(i)/20);
    betah1(i)  =  1/(exp((30-V1(i))/10) + 1);
    alphan1(i) = 0.01*(10-V1(i))/(exp((10-V1(i))/10) - 1);
    betan1(i)  =  0.125*exp(-V1(i)/80);
        
    m1(i) = alpham1(i)/(alpham1(i) + betam1(i));
    h1(i) = alphah1(i)/(alphah1(i) + betah1(i));
    n1(i) = alphan1(i)/(alphan1(i) + betan1(i));


    V2(i) = (-R/2)*((Iext - ((2*V1(i))/R) -((GNa*(m1(i)^3)*h1(i)*(V1(i)-VNa)) + (GK*(n1(i)^4)*(V1(i)-VK))+(GL*(V1(i)-VL))))/c); 


    alpham2(i) = (0.1*(25-V2(i)))/(exp((25-V2(i))/10) - 1);
    betam2(i) =  4*exp(-V2(i)/18);
    alphah2(i) = 0.07*exp(-V2(i)/20);
    betah2(i)  =  1/(exp((30-V2(i))/10) + 1);
    alphan2(i) = 0.01*(10-V2(i))/(exp((10-V2(i))/10) - 1);
    betan2(i)  =  0.125*exp(-V2(i)/80);
        
    m2(i) = alpham2(i)/(alpham2(i) + betam2(i));
    h2(i) = alphah2(i)/(alphah2(i) + betah2(i));
    n2(i) = alphan2(i)/(alphan2(i) + betan2(i));


    V4(i) = -(R)*((((-2*V2(i) + V1(i))/R) -((GNa*(m2(i)^3)*h2(i)*(V2(i)-VNa)) + (GK*(n2(i)^4)*(V2(i)-VK))+(GL*(V2(i)-VL))))/c); 


    alpham4(i) = (0.1*(25-V4(i)))/(exp((25-V4(i))/10) - 1);
    betam4(i) =  4*exp(-V4(i)/18);
    alphah4(i) = 0.07*exp(-V4(i)/20);
    betah4(i)  =  1/(exp((30-V4(i))/10) + 1);
    alphan4(i) = 0.01*(10-V4(i))/(exp((10-V4(i))/10) - 1);
    betan4(i)  =  0.125*exp(-V4(i)/80);
        
    m4(i) = alpham4(i)/(alpham4(i) + betam4(i));
    h4(i) = alphah4(i)/(alphah4(i) + betah4(i));
    n4(i) = alphan4(i)/(alphan4(i) + betan4(i));


    V6(i) = -R*((((-2*V4(i) + V2(i))/R) -((GNa*(m4(i)^3)*h4(i)*(V4(i)-VNa)) + (GK*(n4(i)^4)*(V4(i)-VK))+(GL*(V4(i)-VL))))/c); 
    
    
    alpham6(i) = (0.1*(25-V6(i)))/(exp((25-V6(i))/10) - 1);
    betam6(i) =  4*exp(-V6(i)/18);
    alphah6(i) = 0.07*exp(-V6(i)/20);
    betah6(i)  =  1/(exp((30-V6(i))/10) + 1);
    alphan6(i) = 0.01*(10-V6(i))/(exp((10-V6(i))/10) - 1);
    betan6(i)  =  0.125*exp(-V6(i)/80);
        
    m6(i) = alpham6(i)/(alpham6(i) + betam6(i));
    h6(i) = alphah6(i)/(alphah6(i) + betah6(i));
    n6(i) = alphan6(i)/(alphan6(i) + betan6(i));


    Vfour(i) = R*(-1/2)*((((-2*V6(i))/R) -((GNa*(m6(i)^3)*h6(i)*(V6(i)-VNa)) + (GK*(n6(i)^4)*(V6(i)-VK))+(GL*(V6(i)-VL))))/c); 
     
end

figure,plot(V1,Vfour,V1,V4)