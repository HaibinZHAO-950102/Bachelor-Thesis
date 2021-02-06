clc
clear

N=100;                      %�ڵ���
P0=zeros(1,N+1);            %��ʼѹ��P0����
P=zeros(1,N+1);             %ѹ��P����
DP=zeros(1,N+1);            %��P
XP=zeros(1,N+1);            %�����ڵ�ֵX
H=zeros(1,N+1);             %Ĥ��H
F=zeros(1,N+1);             %F����
B=zeros(1,N+1);             %B����
DF=zeros(1,N+1);            %��F
DL=zeros(1,N+1);            %��CIJʱ��һ���м����
MIDU=zeros(1,N+1);          %�ܶ�ѹ����ϵ
NIANDU=zeros(1,N+1);        %ճ��ѹ����ϵ
CIJ=zeros(N+1,N+1);         %C����
KK=zeros(N,N);              %���յ�K����N*N
KK1=zeros(N+1,N+1);         %���յ�K1����
KK2=zeros(N+1,N+1);         %���յ�K2����
LNODS=zeros(N+1,2);         %�м����

%parameter
w=50;
u=0.2;
NIANDU0=0.08;                    %��ʼճ��
R=0.2;             %�����뾶
G=5.076923076923077e+03;         %���ϲ����������ٻ���
A0=2.2e-8;                       %���ϲ���
E0=G/A0;                         %��������ģ��
U=u*NIANDU0/E0/R;                %�ٶȲ���
W=w/E0/R^2;                      %��λ�����غɣ������ٻ���
B=(sqrt(8.0*W/pi))*R;            %�Ӵ������
PH=E0*B/(4*R);                   %���Ӵ�Ӧ��
T=3.0*(pi^2)*U/(4*W^2);          %��ŵ��������ϵ��
H0=0;                            %��ʼĤ��
XP(1)=-2.5;
XP(N+1)=1;
D=(XP(N+1)-XP(1))/N;             %����
for i= 2:N
    XP(i)=XP(1)+(i-1)*D;         %���ڵ�ֵ
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
        film_thick;          %��Ĥ�񷽳�
        matrix_final;        %������վ���
       DP(2:N+1)=KK\DF(2:N+1)'; %��ⷽ��
       DH0=DP(N+1);
       DP(N+1)=0;
       sum2=0;
       sum3=0;
       for i=1:N
           sum2=sum2+abs(DP(i));
           sum3=sum3+P(i);
       end
       EP=sum2/sum3; 
       P=P+0.4*DP;    %��ɽ������֤�����ľ���ֵ�ȶ��½����ӿ������ٶ�
   end
    H0=H0+0.3*DH0;
end
%dangliangthick;        %��HOIL

figure(1);
plot(H,'r-');    %����Ĥ����״����
hold on;
plot(P,'b-');    %����ѹ���ֲ�����
hold on;
title('��Ĥ��״��ѹ���ֲ�');
axis([1 N 0 1.5]);
set(gca,'Xtick',[1:5:100],'Ytick',[0:0.2:5]);     %����X��������0��100�����5;
                                                    %Y�������0��5�����0.2

grid on;                                            %������ָ���
xlabel('�������');
ylabel('��Ĥ��/ѹ��ֵ');
legend('��Ĥ��H','ѹ��ֵP',-1);
hold off;

DATA=[[1:N+1]' H' P'];
xlswrite('oil.xlsx',DATA)
        