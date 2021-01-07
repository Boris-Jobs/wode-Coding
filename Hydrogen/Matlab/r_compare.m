clc
clear

a0=1;%²£¶û°ë¾¶
r=0:a0/5:50*a0;
[m,n1]=size(r);
for n=1:4
    l=0;
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
    subplot(3,1,1),plot(r,R);
    hold on;
end
legend('n=1,l=0','n=2,l=0','n=3,l=0','n=4,l=0')

for n=2:4
    l=1;
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
    subplot(3,1,2),plot(r,R);
    hold on;
end
legend('n=2,l=1','n=3,l=1','n=4,l=1')

for n=3:4
    l=2;
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
    subplot(3,1,3),plot(r,R);
    hold on;
end
n=4;
l=3;
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
subplot(3,1,3),plot(r,R);
hold on;
legend('n=3,l=2','n=4,l=2','n=4,l=3')

