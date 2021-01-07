clc
clear
N=4;
for n=1:N
    for l=0:n-1
            a0=1;%²£¶û°ë¾¶
            r=0:a0/5:50*a0;
            [m,n1]=size(r);
            a=n-l-1;
            b=factorial(n+l);
            A=sqrt(((2/(n*a0))^3)*factorial(a)/(2*n*b));
            B=zeros(m,n1);
            D=2*r/(n*a0);
            for v=0:a
                sumi=((-1)^(v+1))*(b/(factorial(a-v)*factorial(2*l+1+v)*factorial(v))).*(D.^v);
                B=sumi+B;
            end
            C=exp((-1*r)/(n*a0));
            R=(A.*C.*(D.^l).*B.*r).^2;
            subplot(N,N,(n-1)*N+l+1),plot(r,R),title({['n=',num2str(n)],['l=',num2str(l)]});
    end
end