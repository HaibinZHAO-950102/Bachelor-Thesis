clc
clear
%参数
e=2.1*10^11;             %Pa
miu=0.3;
Ela=e/(1-miu^2);
pusai=atan(5/pi/25);
dw=3.175/1000;           %m
rs=1.75/1000;            %m
aerfa=pi/4;
d0=25/1000;              %m
w=15;                    %N
rou1=2/dw;
rou21=-1/rs;
rou22=2*cos(aerfa)*cos(pusai)/(d0-dw*cos(aerfa));
rou2=rou21+rou22;
R=1/(rou1+rou2);
AA=rou1+rou2;
BB=0;
A=(AA-BB)/2;
B=(AA+BB)/2;
a=(3*w/2/Ela*R)^(1/3);
b=(3*w/2/Ela*R)^(1/3);
po=3*w/2/pi/a/b;
viscosity=0.08;         %Pa s
u=0.01;                  %m/s
stepx=6*b/50;
stepy=6*b/50;
x=-3.5*b:stepx:2.5*b;
y=-2.5*a:stepy:2.5*a;
X=x./b;
Y=y./a;
STEPX=stepx/b;
STEPY=stepy/a;


for i=1:length(x)
    for j=1:length(y)
        if (1-X(i)^2-Y(j)^2<0)
            P0(i,j)=0;
        elseif (1-X(i)^2-Y(j)^2>=0)
            P0(i,j)=(1-X(i)^2-Y(j)^2)^0.5;
        end
    end
end
%figure
%surf(X,Y,P0')


for kkk=1:10

for i=1:length(X)
    for j=1:length(Y)
        delta(i,j)=0;
        for m=1:length(X)
            for n=1:length(Y)
                if i==m&&j==n
                    v(m,n)=0;
                elseif i~=m||j~=n
                    v(m,n)=2*po*b/pi/Ela*(P0(m,n)/sqrt((X(i)-X(m))^2+(Y(j)-Y(n))^2))*STEPX*STEPY;
                end
                delta(i,j)=delta(i,j)+v(m,n);
            end
        end
    end
end
DELTA=delta.*R./b^2;
%figure
%surf(X,Y,DELTA')

HC=2.69*(2.2*10^(-8)*Ela)^0.53*(viscosity*u/Ela/R)^0.67/(w/Ela/R^2)^0.067*(1-0.61*exp(-0.73));
for i=1:length(X)
    for j=1:length(Y)
        H(i,j)=HC+X(i)^2/2+Y(j)^2/2+DELTA(i,j);
    end
end
%figure
%surf(X,Y,H')



%计算方程参数
for i=2:length(X)-1
    for j=2:length(Y)-1
        A(i,j)=1/viscosity*b^4/u*po/R^3*H(i,j)^3;
        B(i,j)=A(i,j);
        C(i,j)=1/viscosity*b^4/u*po/R^3*(((H(i+1,j)^3)-(H(i-1,j)^3))/2/STEPX);
        D(i,j)=1/viscosity*b^4/u*po/R^3*(((H(i,j+1)^3)-(H(i,j-1)^3))/2/STEPY);
        E(i,j)=12*b/R*((H(i+1,j))-(H(i-1,j)))/2/STEPX;
    end
end

%建立系数矩阵
for i=2:length(X)-1
    for j=2:length(Y)-1
        K(i,j)=2*(A(i,j)/STEPX^2+B(i,j)/STEPY^2);
        G(i,j)=-E(i,j)/K(i,j);
        CW(i,j)=(A(i,j)/STEPX^2-C(i,j)/2/STEPX)/K(i,j);
        CE(i,j)=(A(i,j)/STEPX^2+C(i,j)/2/STEPX)/K(i,j);
        CS(i,j)=(B(i,j)/STEPY^2-D(i,j)/2/STEPY)/K(i,j);
        CN(i,j)=(B(i,j)/STEPY^2+D(i,j)/2/STEPY)/K(i,j);
        P(i,j)=CN(i,j)*P0(i,j+1)+CS(i,j)*P0(i,j-1)+CE(i,j)*P0(i+1,j)+CW(i,j)*P0(i-1,j)+G(i,j);
    end
end
for i=1:length(X)
    P(i,1)=0;
    P(i,length(Y))=0;
end
for j=2:length(Y)-1
    P(1,j)=0;
    P(length(X),j)=0;
end
for i=1:length(X)
    for j=1:length(Y)
        if P(i,j)<0
            P(i,j)=0;
        end
    end
end
figure
surf(X,Y,P')
figure

if sum(sum(abs(P0-P)))/sum(sum(P))<0.001
    
    
    
    return
end

kkk
P0=P+0.1*(P0-P);
surf(X,Y,P0')
end


