       % manal{1,k}(i+1) = (((Iext-((GK*((nn{1,k}(i))^4))*(VV{1,k}(i)-VK)))-(GL*(VV{1,k}(i))))/(GNa*hh{1,k}(i)*(VV{1,k}(i)-VNa)))^(1/3);
       % nanal{1,k}(i+1) = ((Iext - ((GNa*(mm{1,k}(i)^3))*(hh{1,k}(i))*(VV{1,k}(i)-VNa)) - (GL*(VV{1,k}(i))))/(((GK*(VV{1,k}(i)-VK)))^(1/4)));
       % hanal{1,k}(i+1) = (Iext-((GK*(nn{1,k}(i)^4))*(VV{1,k}(i)-VK))-(GL*(VV{1,k}(i))))/(GNa*(((mm{1,k}(i))^3)*(VV{1,k}(i)-VNa)));
        clear all
        clc
       
        VV = 0:0.01:120;
        
        Iext = 5;
        GNa = 120 ;
        GK = 36;
        GL = 0.3;
        R = 10;
        VNa = 115;
        VK = -12;
        VL = 10.5995;
        c = 1;
    
        alpham  = zeros(1,12001);
        alphah  = zeros(1,12001);
        alphan  = zeros(1,12001);
        betam   = zeros(1,12001);
        betah   = zeros(1,12001);
        betan   = zeros(1,12001);
        m2      = zeros(1,12001);
        h2      = zeros(1,12001);
        n2      = zeros(1,12001);
        m1anal  = zeros(1,12001);
        h1anal  = zeros(1,12001);
        n1anal  = zeros(1,12001);
        
        for i = 1:12001
        
        alpham(i) = (0.1*(25-VV(i)))/(exp((25-VV(i))/10) - 1);
        betam(i) =  4*exp(-VV(i)/18);
        alphah(i) = 0.07*exp(-VV(i)/20);
        betah(i) =  1/(exp((30-VV(i))/10) + 1);
        alphan(i) = 0.01*(10-VV(i))/(exp((10-VV(i))/10) - 1);
        betan(i) =  0.125*exp(-VV(i)/80);
        
        m2(i)  = (alpham(i)/(alpham(i) + betam(i)));
        h2(i)=  (alphah(i)/(alphah(i) + betah(i)));
        n2(i)=  (alphan(i)/(alphan(i) + betan(i)));
        
        n4anal(i) = ((Iext - ((GNa*((m2(i))^3))*(h2(i))*(VV(i)-VNa)) - (GL*(VV(i)- VL)))/(GK*(VV(i)-VK)));
        h1anal(i) = (Iext-(GK*((n2(i))^4)*(VV(i)-VK))-(GL*(VV(i)- VL )))/(GNa*(((m2(i))^3)*(VV(i)-VNa)));
        m3anal(i) = ((Iext-((GK*(n2(i)^4))*(VV(i)-VK)))-(GL*(VV(i)- VL)))/(GNa*h2(i)*(VV(i)-VNa));
        m1anal(i) = nthroot(m3anal(i),3);
        n1anal(i) = (n4anal(i))^(1/4);
        fpm(i)    = m2(i)- m1anal(i);
        fph(i)    = h2(i)- h1anal(i);
        fpn(i)    = n2(i)- n1anal(i);
        end
        %figure,plot(m1anal,VV)
        %figure,plot(h1anal,VV)
        %figure,plot(n1anal,VV)
        %figure,plot(VV,fpm)
        %figure,plot(VV,fph)
        %figure,plot(VV,fpn)
        %figure,plot(m3anal,VV)
        %figure,plot(m2-m1anal)
        figure,plot(VV,m2,VV,m1anal)
        figure,plot(VV,n2,VV,n1anal)
        figure,plot(VV,h2,VV,h1anal)