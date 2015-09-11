function z=Analyzer(self,x)
%Analyzer is a method of the class FPDAnalyzer which calculates the
%approximated probability of given protein level(s).
%  z=obj.Analyzer(x) returns the approximated probability of protein
%  level(s) given in the one dimensional nonnegative integer array x for
%  the model defined in obj. 
    u=unique(x);
    x=x(:)';
    u=u(:)';
    t=diag(0:0.1:1)*ones(11,length(u))+ones(11,length(u))*diag(u);
    s=size(t);
    t=t(:);
    z=self.fun(t);
    z=reshape(z,s);
    z1=z(1:end-1,:);
    z2=z(2:end,:);
    zt=sum(z1+z2,1)*0.05;
    z=x;
    for i=1:length(u)
        if(~isreal(zt(i))||isinf(zt(i)))
            zt(i)=integral(@(v)self.fun(v),u(i),u(i)+1);
        end
        z(x==u(i))=zt(i);
    end
end
