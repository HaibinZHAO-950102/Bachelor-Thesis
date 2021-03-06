function x=gauss(a,b)  %用途：顺序Gauss消去法解线性方程组ax=b
n=length(b);
a=[a,b];
for k=1:(n-1)   %消元过程
 a((k+1):n,(k+1):(n+1))=a((k+1):n,(k+1):(n+1))-a((k+1):n,k)/a(k,k)*a(k,(k+1):(n+1));
 a((k+1):n,k)=zeros(n-k,1);
end
x=zeros(n,1);
x(n)=a(n,n+1)/a(n,n);
for k=n-1:-1:1
    x(k,:)=a(k,(n+1))-a(k,(k+1):n)*x((k+1):n)/a(k,k);
end