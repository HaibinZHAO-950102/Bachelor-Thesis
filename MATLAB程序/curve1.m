clc
clear
%计算反向器XYZ坐标
m=4/3*pi*(3.175/2/10)^3*7.85/1000;
n=1/350;
rs=100;
pusai=(39/60+3)/180*pi;
P=5;
R=12.5;
beita=40/180*pi;
h=5-beita/2/pi*5;
a=2.2;
sita1=-beita/2:n:beita/2;
A=2*pi*h/(2*pi*h+beita*P);
B=h/A/beita;
x=(-beita/2/pi*sin((sita1+beita/2)*2*pi/beita)+A*(sita1+beita/2))*B;
x=x;
l=a/2*((1-cos(2*pi*(sita1+beita/2)/beita)-1/4*(1-cos(4*pi*(sita1+beita/2)/beita))));
r=R+l;
y=r.*sin(sita1);
z=r.*cos(sita1);
plot3(z,y,x)

%计算螺旋线坐标
sita2=beita/2:n:2*pi-beita/2;
X=sita2/2/pi*5;
X=fliplr(X);
X=X-beita/4/pi*5;
Y=R.*sin(sita2);
Z=R.*cos(sita2);
hold on
plot3(Z,Y,X);xlabel('Z');ylabel('Y');zlabel('X');axis equal;
x=[x X];y=[y Y]; z=[z Z];
curve=[x' y' z'];
xlswrite('curve.xlsx',curve)