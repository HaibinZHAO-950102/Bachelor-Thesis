clc
clear
%∏’∂»º∆À„
N=10;
for Q=0:N:10000
q=Q/66;
dw=3.175;
d0=25;
rs=1.75;
alfa=pi/4;
pusai=atan(5/pi/25);
e=2.1*10^5;
miu=0.3;
rou11=2/dw;
rou12=2/dw;
rou21=-1/rs;
rou22=-2*cos(alfa)*cos(pusai)/(d0-dw*cos(alfa));
E=e/(1-miu^2);
rou1=rou11+rou22;
rou2=rou12+rou21;
rou=rou1+rou2;
k=1.0339*(rou1/rou2)^0.636;
sigema=1.0003+0.5968*(rou1/rou2);
T=1.5277+0.6023*(log((rou1/rou2))/log(exp(1)));
ax=(2*k^2*sigema/pi)^(1/3);
bx=(2*sigema/pi/k)^(1/3);
dx=2*T/pi*(pi/2/k^2/sigema)^(1/3);
a=ax*(3*q/E/rou)^(1/3);
b=bx*(3*q/E/rou)^(1/3);
d=dx*(3*q/E/rou)^(2/3)*rou/2;
roun=rou11+rou12+rou21+rou22;
rous=rou11+rou12+rou21-rou22;
dn=1/pi*(3*q/E)^(2/3)*0.209940957*(roun^(1/3)+rous^(1/3));
da(Q/N+1)=dn*cos(pusai)/sin(alfa);
K(Q/N+1)=Q/da(Q/N+1)/10^6;
end
F=0:N:Q;
K(1)=0;
da(1)=0;
plot(F,K)
xlswrite('stiffness.xlsx',[F' K'])


da=da(1:100);
F=F(1:100);
n1=find(F==200);
n2=find(F==500);
dan1=-1.*da+2*da(n1);
dan2=-1.*da+2*da(n2);
figure
plot(da,F,dan1,F,dan2,F)
xlswrite('preload.xlsx',[da' F' dan1' F' dan2' F'])

