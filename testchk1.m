clear all
clc
[x,y] = meshgrid ([-4:0.01:4]);
a = 2;
Z = exp((-(((sqrt(x.*x + y.*y)))- 1).^2)/(a^2));
figure,mesh(Z)