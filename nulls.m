clear all
clc
[x,y] = meshgrid ([-4:0.005:4]);
a = 20;
Z = exp(-(((sqrt(x.*x + y.*y))-1)^2)/(2*(a^2)));
figure,mesh(Z)