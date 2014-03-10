% Author: Housam Binous

% The Cubic Map

% National Institute of Applied Sciences and Technology, Tunis, TUNISIA

% Email: binoushousam@yahoo.com

% Reference :
% Varma and Morbidelli, Mathematical Methods in Chemical Engineering, 
% Oxford University Press, 1997.

global ra

% plotting the bifurcation diagram

r=0:0.1:8;

figure(1)
hold on
for i=1:81;
    ra=r(i);
    a=0.1;
    for j=1:900
    a=feval(@cubic2,a);
    end;
    a(1)=a;
    for j=2:100
    a(j)=feval(@cubic2,a(j-1));
    end
    rp=ra*ones(1,100);
    plot(rp,a,'r.')
end;

% A cycle of Period 3 appears at r?6.05

figure(2)
ra=6.05;
a=0.2;
for j=1:200
    a=feval(@cubic2,a);
    end;
    a(1)=a;
    for j=2:100
    a(j)=feval(@cubic2,a(j-1));
    end;
plot(a,'r.')

% Chaotic behavior appears at r>6.05

figure(3)
ra=6.8;
a=0.2;
for j=1:200
    a=feval(@cubic2,a);
    end;
    a(1)=a;
    for j=2:100
    a(j)=feval(@cubic2,a(j-1));
    end;
plot(a,'b.')
axis([0 100 0 1.2])
