  for i=1:N+1                                                      %            
            H(i)=H0+XP(i)^2/2;        %Ĥ��ʽǰ����                  %
            DH=0;                                                        %
            for j=2:N+1                                                  %
                DL(j-1)=XP(j)-XP(j-1);                                   %
                if(i==j||i==j-1)
                    aa=XP(i);
                    bb=XP(j)+0.1*DL(j-1);
                    if aa==bb
                        CIJ(i,j)=(-0.5/3.14)*DL(j-1)*log(0.01);
                    else            
                        CIJ(i,j)=(-0.5/3.14)*DL(j-1)*log((aa-bb)^2);          
                    end
                else                                                     %
                    CIJ(i,j)=(-0.5/3.14)*DL(j-1)*log((XP(i)-XP(j))^2);     %                    
                end                                                      %
                DH=DH+CIJ(i,j)*P(j);    %Ĥ��ʽ�󲿷֦�(X)             %
            end                                                          %
            H(i)=H(i)+DH;               %Ĥ�񡡡�                        %
        end                                                              %  