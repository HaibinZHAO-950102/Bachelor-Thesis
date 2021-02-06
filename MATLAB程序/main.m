clc
clear

N=100;                      %节点数
P0=zeros(1,N+1);            %初始压力P0向量
P=zeros(1,N+1);             %压力P向量
DP=zeros(1,N+1);            %ΔP
XP=zeros(1,N+1);            %各个节点值X
H=zeros(1,N+1);             %膜厚H
F=zeros(1,N+1);             %F向量
B=zeros(1,N+1);             %B向量
DF=zeros(1,N+1);            %ΔF
DL=zeros(1,N+1);            %求CIJ时的一个中间变量
MIDU=zeros(1,N+1);          %密度压力关系
NIANDU=zeros(1,N+1);        %粘度压力关系
CIJ=zeros(N+1,N+1);         %C矩阵
KK=zeros(N,N);              %最终的K矩阵N*N
KK1=zeros(N+1,N+1);         %最终的K1矩阵
KK2=zeros(N+1,N+1);         %最终的K2矩阵
LNODS=zeros(N+1,2);         %中间矩阵

%parameter
w=50;
u=0.2;
NIANDU0=0.08;                    %初始粘度
R=0.2;             %当量半径
G=5.076923076923077e+03;         %材料参数（无量纲化）
A0=2.2e-8;                       %材料参数
E0=G/A0;                         %当量弹性模量
U=u*NIANDU0/E0/R;                %速度参数
W=w/E0/R^2;                      %单位长度载荷（无量纲化）
B=(sqrt(8.0*W/pi))*R;            %接触区半宽
PH=E0*B/(4*R);                   %最大接触应力
T=3.0*(pi^2)*U/(4*W^2);          %雷诺方程右项系数
H0=0;                            %初始膜厚
XP(1)=-2.5;
XP(N+1)=1;
D=(XP(N+1)-XP(1))/N;             %步长
for i= 2:N
    XP(i)=XP(1)+(i-1)*D;         %各节点值
end   


for i=1:N+1
    if(XP(i)^2<=1)
        P0(i)=sqrt(1-XP(i)^2);
    end
end
P=P0;
DH0=1;
MAXIT=100;
ITER=0;
while (abs(DH0)>0.01&&ITER<=MAXIT)
    ITER=ITER+1;
    IT=0;
    EP=1;
    while (EP>0.01&&IT<=MAXIT)
        IT=IT+1;
        film_thick;          %求膜厚方程
        matrix_final;        %求解最终矩阵
       DP(2:N+1)=KK\DF(2:N+1)'; %求解方程
       DH0=DP(N+1);
       DP(N+1)=0;
       sum2=0;
       sum3=0;
       for i=1:N
           sum2=sum2+abs(DP(i));
           sum3=sum3+P(i);
       end
       EP=sum2/sum3; 
       P=P+0.4*DP;    %下山法，保证函数的绝对值稳定下降，加快收敛速度
   end
    H0=H0+0.3*DH0;
end
%dangliangthick;        %求HOIL

figure(1);
plot(H,'r-');    %划出膜厚形状曲线
hold on;
plot(P,'b-');    %划出压力分布曲线
hold on;
title('润滑膜形状和压力分布');
axis([1 N 0 1.5]);
set(gca,'Xtick',[1:5:100],'Ytick',[0:0.2:5]);     %设置X轴的坐标从0到100，间距5;
                                                    %Y轴坐标从0到5，间距0.2

grid on;                                            %画坐标分隔线
xlabel('Ｘ坐标点');
ylabel('润滑膜厚/压力值');
legend('润滑膜厚H','压力值P',-1);
hold off;

DATA=[[1:N+1]' H' P'];
xlswrite('oil.xlsx',DATA)
        