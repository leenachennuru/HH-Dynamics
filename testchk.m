clear all
clc
[x,y] = meshgrid (-4:0.01:4);
a = 20;
Z = (x.^2 + y.^2)-1;
figure,mesh(Z)