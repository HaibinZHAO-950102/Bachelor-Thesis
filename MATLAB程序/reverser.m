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
sita1=0:n:beita;
A=2*pi*h/(2*pi*h+beita*P);
B=h/A/beita;
x=(-beita/2/pi*sin((sita1)*2*pi/beita)+A*(sita1))*B;
x=x+P-h;
l=a/2*((1-cos(2*pi*(sita1)/beita)-1/4*(1-cos(4*pi*(sita1)/beita))));
r=R+l;
y=r.*sin(sita1);
z=r.*cos(sita1);
%plot3(z,y,x)

%计算螺旋线坐标
sita2=beita:n:2*pi;
X=sita2/2/pi*5;
X=fliplr(X);
Y=R.*sin(sita2);
Z=R.*cos(sita2);
%hold on
%plot3(Z,Y,X);xlabel('Z');ylabel('Y');zlabel('X');axis equal;


%计算路径长度
L=0;
curve=curve./1000;
x=x./1000;
y=y./1000;
z=z./1000;
for i=1:length(x)-1
   L=L+sqrt((curve(i+1,1)-curve(i,1))^2+(curve(i+1,2)-curve(i,2))^2+(curve(i+1,3)-curve(i,3))^2);
end
L=L+sqrt((curve(1,1)-curve(length(x)-1,1))^2+(curve(1,2)-curve(length(x)-1,2))^2+(curve(1,3)-curve(length(x)-1,3))^2);
C=1/(13.6225/2/12.5);

sita=[sita1 sita2];

%计算XYZ方向速度
for i=1:length(sita1)-2
    vx(i)=(x(i+1)-x(i))/((sita(i+1)-sita(i))/(2*pi)/rs);
end
for i=length(sita1)-1:length(sita)
    vx(i)=-5/1000*rs;
end
%figure;plot(sita,vx)
VX=[sita'*C vx'/C];
xlswrite('vx.xlsx',VX);

for i=1:length(sita1)
    vy(i)=(y(i+1)-y(i))/((sita(i+1)-sita(i))/(2*pi)/rs);
end
for i=length(sita1)+1:length(sita)
    vy(i)=R/1000*rs*2*pi*cos(sita(i));
end
%figure;plot(sita,vy);
VY=[sita'*C vy'/C];
xlswrite('vy.xlsx',VY);

for i=1:length(sita1)
    vz(i)=(z(i+1)-z(i))/((sita(i+1)-sita(i))/(2*pi)/rs);
end
for i=length(sita1)+1:length(sita)
    vz(i)=-R/1000*rs*2*pi*sin(sita(i));
end
%figure;plot(sita,vz)
VZ=[sita'*C vz'/C];
xlswrite('vz.xlsx',VZ);

%计算XYZ方向加速度
for i=1:length(sita1)-2
    ax(i)=(vx(i+1)-vx(i))/((sita(i+1)-sita(i))/(2*pi)/rs);
end
for i=length(sita1)-1:length(sita)
    ax(i)=0;
end
%figure;plot(sita,ax)
Ax=[sita'*C ax'/C];
xlswrite('ax.xlsx',Ax);

for i=1:length(sita1)-2
    ay(i)=(vy(i+1)-vy(i))/((sita(i+1)-sita(i))/(2*pi)/rs);
end
for i=length(sita1)-1:length(sita)
    ay(i)=-R/1000*rs*2*pi*rs*2*pi*sin(sita(i));
end
%figure;plot(sita,ay)
Ay=[sita'*C ay'/C];
xlswrite('ay.xlsx',Ay);

for i=1:length(sita1)-2
    az(i)=(vz(i+1)-vz(i))/((sita(i+1)-sita(i))/(2*pi)/rs);
end
for i=length(sita1)-1:length(sita)
    az(i)=-R/1000*rs*2*pi*rs*2*pi*cos(sita(i));
end
%figure;plot(sita,az)
Az=[sita'*C az'/C];
xlswrite('az.xlsx',Az);

k=length(sita)/25;
%叠加25个滚珠的加速度
for i=2:25
    for j=k*(i-1)+1:length(x)
        ax(i,j)=ax(1,j-(i-1)*k);
    end
end
for i=2:25
    for j=1:k*(i-1)
        ax(i,j)=ax(1,length(sita)-k*(i-1)+j);
    end
end

for i=2:25
    for j=k*(i-1)+1:length(x)
        ay(i,j)=ay(1,j-(i-1)*k);
    end
end
for i=2:25
    for j=1:k*(i-1)
        ay(i,j)=ay(1,length(sita)-k*(i-1)+j);
    end
end

for i=2:25
    for j=k*(i-1)+1:length(x)
        az(i,j)=az(1,j-(i-1)*k);
    end
end
for i=2:25
    for j=1:k*(i-1)
        az(i,j)=az(1,length(sita)-k*(i-1)+j);
    end
end

for i=1:length(x)
    AX(i)=sum(ax(:,i));
    AY(i)=sum(ay(:,i));
    AZ(i)=sum(az(:,i));
end
%figure;plot(sita,AX);figure;plot(sita,AY);figure;plot(sita,AZ)
FX=[sita'*C AX'/C*m];
FY=[sita'*C AY'/C*m];
FZ=[sita'*C AZ'/C*m];
xlswrite('fx.xlsx',FX);
xlswrite('fy.xlsx',FY);
xlswrite('fz.xlsx',FZ);
