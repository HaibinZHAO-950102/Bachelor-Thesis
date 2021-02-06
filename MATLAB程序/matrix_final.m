 NNODZ=2;                                                                              %
       NIANDU0=0.08;               % ��ʼճ�Ȧ�0                                             %
       for i=1:N+1                 %���ܶȺ�ճ�ȡ���������������                             %
           MIDU(i)=1.0+((0.6e-9)*P(i)*PH)/(1.0+(1.7e-9)*P(i)*PH);                            %
           NIANDU(i)=exp((log(NIANDU0)+9.67)*(-1+(1+(5.1e-9)*P(i)*PH)^0.6));                %
       end                                                                                   %
       KK1=zeros(N+1,N+1);              %��0����                                             %
       KK2=zeros(N+1,N+1);                                                                   %
       F=zeros(1,N+1);                                                                       %
       B=zeros(1,N+1);                                                                       %
       for NEL=2:N+1                                                                         %
           matrix_element;          %������еĸ���Ԫ��
           for i=1:NNODZ                                                                     %
               LNODS(NEL,i)=NEL+i-2;                                                         %
           end                                                                               %
           for i=1:NNODZ                                                                     %
               ISTRST=LNODS(NEL,i);                                                          %
               IELEMT=i;                                                                     %
               for j=1:NNODZ                                                                 %
                   JSTRST=LNODS(NEL,j);                                                      %
                   JELEMT=j;                %��Ӻ�ģ�1                                     %
                   KK1(ISTRST,JSTRST)=KK1(ISTRST,JSTRST)+A(IELEMT,JELEMT);                   %
               end                                                                           %
               F(ISTRST)=F(ISTRST)+FF(IELEMT);     %��Ӻ�ģ�����                           %
               B(ISTRST)=B(ISTRST)+BB(IELEMT);     %��Ӻ�ģ�������                         %
               for j=1:N                                                                     %
                   KK2(ISTRST,j)=KK2(ISTRST,j)+AL(IELEMT,j);  %��Ӻ�ģ�2����               %
               end                                                                           %
           end                                                                               %
       end                                                                                   %
       sum1=0;                                                                               %
       for i=2:N                                                                             %
           sum=0;                                                                            %
           for j=2:N   %�˴������е�N+1*N+1�ľ���ֱ�Ӱ�N*N���㣬�Ա��ڽⷽ�̾���      ��     %
               KK(i-1,j-1)=KK1(i,j)+KK2(i,j);    %��1���2��Ӻ�ģ˾���                     %
               sum=sum+KK(i-1,j-1)*P(j);                                                     %
           end                                                                               %
           DF(i)=F(i)-sum;                    %���ƣ���-��*��0                               %
           KK(i-1,N)=B(i);                    %�������ӵ��˾����С���                        %
           KK(N,i-1)=(XP(i+1)-XP(i-1))/2;     %���������뵽�˾�����                          %
           sum1=sum1+P(i)*(XP(i+1)-XP(i));                                                   %
       end                                                                                   %
       DF(N+1)=pi/2-sum1;                      %����                                         %
       KK(N,N)=0;                                                                            %
       %....................................................                 . ..............%