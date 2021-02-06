 NNODZ=2;                                                                              %
       NIANDU0=0.08;               % 初始粘度η0                                             %
       for i=1:N+1                 %求密度和粘度　　　　　　　　                             %
           MIDU(i)=1.0+((0.6e-9)*P(i)*PH)/(1.0+(1.7e-9)*P(i)*PH);                            %
           NIANDU(i)=exp((log(NIANDU0)+9.67)*(-1+(1+(5.1e-9)*P(i)*PH)^0.6));                %
       end                                                                                   %
       KK1=zeros(N+1,N+1);              %置0　　                                             %
       KK2=zeros(N+1,N+1);                                                                   %
       F=zeros(1,N+1);                                                                       %
       B=zeros(1,N+1);                                                                       %
       for NEL=2:N+1                                                                         %
           matrix_element;          %求矩阵中的各个元素
           for i=1:NNODZ                                                                     %
               LNODS(NEL,i)=NEL+i-2;                                                         %
           end                                                                               %
           for i=1:NNODZ                                                                     %
               ISTRST=LNODS(NEL,i);                                                          %
               IELEMT=i;                                                                     %
               for j=1:NNODZ                                                                 %
                   JSTRST=LNODS(NEL,j);                                                      %
                   JELEMT=j;                %相加后的Ｋ1                                     %
                   KK1(ISTRST,JSTRST)=KK1(ISTRST,JSTRST)+A(IELEMT,JELEMT);                   %
               end                                                                           %
               F(ISTRST)=F(ISTRST)+FF(IELEMT);     %相加后的Ｆ向量                           %
               B(ISTRST)=B(ISTRST)+BB(IELEMT);     %相加后的Ｂ向量　                         %
               for j=1:N                                                                     %
                   KK2(ISTRST,j)=KK2(ISTRST,j)+AL(IELEMT,j);  %相加后的Ｋ2　　               %
               end                                                                           %
           end                                                                               %
       end                                                                                   %
       sum1=0;                                                                               %
       for i=2:N                                                                             %
           sum=0;                                                                            %
           for j=2:N   %此处把书中的N+1*N+1的矩阵直接按N*N计算，以便于解方程矩阵      　     %
               KK(i-1,j-1)=KK1(i,j)+KK2(i,j);    %Ｋ1与Ｋ2相加后的Ｋ矩阵                     %
               sum=sum+KK(i-1,j-1)*P(j);                                                     %
           end                                                                               %
           DF(i)=F(i)-sum;                    %ΔＦ＝Ｆ-Ｋ*Ｐ0                               %
           KK(i-1,N)=B(i);                    %Ｂ向量加到Ｋ矩阵中　　                        %
           KK(N,i-1)=(XP(i+1)-XP(i-1))/2;     %Ｄ向量加入到Ｋ矩阵中                          %
           sum1=sum1+P(i)*(XP(i+1)-XP(i));                                                   %
       end                                                                                   %
       DF(N+1)=pi/2-sum1;                      %ΔＷ                                         %
       KK(N,N)=0;                                                                            %
       %....................................................                 . ..............%