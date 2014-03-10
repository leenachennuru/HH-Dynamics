x = linspace(-3,3,30);
y = linspace(-3,3,30);
Z = cell(30,30);
Z1 = cell(30,30);
for i = 1:30
for j = 1:30
Z{i,j} = [x(i),y(j)];
end
end
for m = 1:30
for n = 1:30
Z1{m,n} = [y(n), (((1/3)-(x(m)^2))*y(n))-((x(m))^5)];
end
end
