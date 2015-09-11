function z=LaguerreApproximate(self,x)
    e=self.ecount;
    e=e(:)';
    ntotal=length(e)-1;
    c=(e(3)-e(2)^2)/e(2);
    a=e(2)/c;
    coef=gso(ntotal,a);
    x=x./c;
    e=e./(c.^(0:ntotal));
    e=e';
    xs=coef*e;
    xx=ones(ntotal+1,length(x))*diag(x);
    xp=diag(0:ntotal)*ones(ntotal+1,length(x));
    xx=xx.^xp;
    xx=coef*xx;
    xs=xs(:)';
    z=xs*xx;
    xu=(x.^(a-1)).*exp(-x)./gamma(a);
    xu=xu(:)';
    z=z.*xu;
    z=z./c;
    z=(z>=0).*z;
end

function z=gso(n,a)
    z=zeros(n+1);
    z(1,1)=1;
    for i=1:n
        xn=zeros(1,n+1);
        xn(i+1)=1;
        yn=xn;
        for j=0:(i-1)
            yn=yn-ip(xn,z(j+1,:),a)*z(j+1,:);
        end
        yn=yn./sqrt(ip(yn,yn,a));
        z(i+1,:)=yn;
    end
end

function z=ip(x,y,a)
    x=x(:)';
    y=y(:)';
    z=conv(x,y);
    w=gamma(a+(0:(length(z)-1)))./gamma(a);
    w=w(:)';
    z=sum(z.*w);
end
