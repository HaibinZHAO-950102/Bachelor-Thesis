           A=zeros(2,2);     %K1                                                 %           %
           AL=zeros(2,N+1);  %K2                                                 %           %
           BB=linspace(0,0,2);                                                   %           %
           FF=linspace(0,0,2);                                                   %           %
           E=linspace(0,0,N+1);                                                  %           %
           for i=NEL-1:NEL                                                       %           %
               E(i)=MIDU(i)*H(i)^3/NIANDU(i);                                    %           %
           end                                                                   %           %
           for i=1:N              %K2                                            %           %
               AL(1,i)=T*(MIDU(NEL-1)*CIJ(NEL-1,i)+MIDU(NEL)*CIJ(NEL,i))/2;      %           %
               AL(2,i)=-AL(1,i);                                                 %           %
           end                                                                   %           %
           A(1,1)=(E(NEL-1)+E(NEL))/(2.0*D);        %K1                            %           %
           A(1,2)=-A(1,1);                                                       %           %
           A(2,1)=A(1,2);                                                        %           %
           A(2,2)=A(1,1);                                                        %           %
           FF(1)=-T*(MIDU(NEL-1)+MIDU(NEL))/2.0*(H0+(XP(NEL)^3-XP(NEL-1)^3)/(6*D)); %          %         
           FF(2)=-FF(1);                                                         %           %
           for i=NEL-1:NEL                                                       %           %
               if(i==1||i==N+1)                                                  %           %
                   DP(i)=0;                                                      %           %
               else                                                              %           %
                   DP(i)=(P(i+1)-P(i-1))/(XP(i+1)-XP(i-1));                      %           %
               end                                                               %           %
           end                                                                   %           %
           BB(1)=-3.0*(MIDU(NEL-1)*H(NEL-1)^2/NIANDU(NEL-1)*DP(NEL-1)+MIDU(NEL)*... %          %
               H(NEL)^2/NIANDU(NEL)*DP(NEL))/2+T*(MIDU(NEL-1)+MIDU(NEL))/2.0;      %           %
           BB(2)=-BB(1);                                          