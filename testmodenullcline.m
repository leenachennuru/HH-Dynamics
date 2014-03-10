clear all
h = 0.596; 
n = 0.318;
GNa = 120 ;
GK = 36;
GL=0.3;
R = 10;
VNa = 115;
VK = -12;
VL=10.5995;
c = 1;
Iext = 0.20;
A = zeros(1,11401);
y1 = zeros(1,11401);
B = 1:0.01:115;
for i = 1:11401
A(i) = (((20-((36*(0.318^4))*(B(i)+12)))-(0.3*(B(i)-10.5995)))/(120*0.596*(B(i)-115)))^(1/3);
a(i) = (0.1*(25-B(i)))/(exp((25-B(i))/10) - 1);
b(i) =  4*exp(-B(i)/18);
y1(i) = a(i)/(a(i)+ b(i));
end
figure,plot(B,y1,B,A)
figure,plot(B,y1)
figure,plot(B,A)

