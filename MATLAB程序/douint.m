clc
clear
step=0.002;
x=-1:step:1;
y=-1:step:1;
for i=1:length(x)
    for j=1:length(y)
        z(i,j)=2-(x(i).^2+y(j).^2);
    end
end
surf(x,y,z)
V=0;
for i=1:length(x)-1
    for j=1:length(y)-1
        V=V+z(i+1,j+1)*step*step;
    end
end
