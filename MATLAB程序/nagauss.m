function x=gauss(a,b)  %��;��˳��Gauss��ȥ�������Է�����ax=b
n=length(b);
a=[a,b];
for k=1:(n-1)   %��Ԫ����
 a((k+1):n,(k+1):(n+1))=a((k+1):n,(k+1):(n+1))-a((k+1):n,k)/a(k,k)*a(k,(k+1):(n+1));
 a((k+1):n,k)=zeros(n-k,1);
end
x=zeros(n,1);
x(n)=a(n,n+1)/a(n,n);
for k=n-1:-1:1
    x(k,:)=a(k,(n+1))-a(k,(k+1):n)*x((k+1):n)/a(k,k);
end