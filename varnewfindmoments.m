function z=varnewfindmoments(self)
    para=self.para;
    ntotal=self.order;
    dtotal=1;
    t0=para(1);t1=para(2);k1=para(3);kp=para(4);dm=para(5);dp=para(6);
    e=zeros(dtotal+1,ntotal+1,ntotal+1);
    e(:,1,1)=getinit(dtotal,t0,t1);
    for n=1:ntotal
        e=recuserve(t0,t1,k1,kp,dm,dp,e,dtotal,n);
    end
    z=zeros(1,ntotal+1);
    for i=0:ntotal
        z(i+1)=sum(e(:,1,i+1));
    end
end

function z=getinit(d,t0,t1)
    c=zeros(d+1);c(end,:)=1;b=zeros(d+1,1);b(end)=1;
    delta=0;c(1,delta+1)=-(t0*delta+t1*(d-delta));c(1,delta+2)=t0*(delta+1);
    for delta=1:(d-1)
        c(delta+1,delta)=t1*(d-delta+1);
        c(delta+1,delta+1)=-(t0*delta+t1*(d-delta));
        c(delta+1,delta+2)=t0*(delta+1);
    end
    z=c\b;
end

function z=recuserve(t0,t1,k1,kp,dm,dp,e,dtotal,n)
    l=zeros((dtotal+1)*(n+1));
    r=zeros(1,(dtotal+1)*(n+1))';
    for d=0:dtotal
        for x=0:n
            y=n-x;
            if(d>0)
                l(d*(n+1)+x+1,(d-1)*(n+1)+x+1)=l(d*(n+1)+x+1,(d-1)*(n+1)+x+1)+t1*(dtotal-d+1);
            end
            if(d<dtotal)
                l(d*(n+1)+x+1,(d+1)*(n+1)+x+1)=l(d*(n+1)+x+1,(d+1)*(n+1)+x+1)+t0*(d+1);
            end
            for k=1:x
                r(d*(n+1)+x+1)=r(d*(n+1)+x+1)+k1*d*nchoosek(x,k)*e(d+1,x-k+1,y+1);
                if(k==1)
                    l(d*(n+1)+x+1,d*(n+1)+x+1)=l(d*(n+1)+x+1,d*(n+1)+x+1)+dm*nchoosek(x,k)*(-1)^k;
                else
                    r(d*(n+1)+x+1)=r(d*(n+1)+x+1)+dm*nchoosek(x,k)*(-1)^k*e(d+1,x-k+2,y+1);
                end
            end
            for k=1:y
                if(k==1)
                    l(d*(n+1)+x+1,d*(n+1)+x+2)=l(d*(n+1)+x+1,d*(n+1)+x+2)+kp*nchoosek(y,k);
                    l(d*(n+1)+x+1,d*(n+1)+x+1)=l(d*(n+1)+x+1,d*(n+1)+x+1)+dp*nchoosek(y,k)*(-1)^k;
                else
                    r(d*(n+1)+x+1)=r(d*(n+1)+x+1)+kp*nchoosek(y,k)*e(d+1,x+2,y-k+1)...
                        +dp*nchoosek(y,k)*(-1)^k*e(d+1,x+1,y+2-k);
                end
            end
            l(d*(n+1)+x+1,d*(n+1)+x+1)=l(d*(n+1)+x+1,d*(n+1)+x+1)-(t0*d+t1*(dtotal-d));
        end
    end
    r=-r;
    z=l\r;
    for d=0:dtotal
        for x=0:n
            y=n-x;
            e(d+1,x+1,y+1)=z(d*(n+1)+x+1);
        end
    end
    z=e;
end
