clc
clear

i=4;
for l=1:i
    for m=1:i
        if m<=l
            theta=-pi/2:0.01:pi/2;
            phi=0:0.02:2*pi;
            [phi2,theta2]=meshgrid(phi,theta);
            a=legendre(l-1,sin(theta2));
            [h,n,o]=size(a);
            s=zeros(n,o);
            for x=1:n
                for y=1:o
                    s(x,y)=a(m,x,y);
                end
            end
            r2=((factorial(l-m)*(2*l-1))/(factorial(l+m)*4*pi)).*(s).^2;
            [x,y,z]=sph2cart(phi2,theta2,r2);
            subplot(i,i,(l-1)*i+m)
            mesh(x,y,z),title({['角量子数：l=',num2str(l-1)],['磁量子数：m=',num2str(m-1)]});
        else
            subplot(i,i,(l-1)*i+m),title({['404 not found'],['41821256 陈哲']});
        end
    end
end



