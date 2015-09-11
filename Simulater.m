function self=Simulater(self,varargin)
%Simulater is a method of the class FPDAnalyzer which using the Gillespie
%algorithm to simulate the reaction path of gene expression from a total
%inactive state (gene inactive, no mRNA and protein).
%  z=obj.Simulater(T,n) simulate a reaction patt of n of the model defined
%  in obj. The program will stop when time T is reached or all the first
%  four moments of the simulated protein population are differ no more than
%  5% from the theoretical ones. the property simudata is updated by the
%  newly generated data.
%  z=obj.Simulater(T) take n=1000 as the default value. 
%  z=obj.Simulater() also take T=1000 as default value.
%  The output stores the updated instance. 
%  The simulation may take a lot of time.

    para1=num2cell(self.para);
    [t0,t1,k1,kp,dm,dp]=deal(para1{:});
    if(isempty(varargin))
        T=Inf;
        n=1000;
    elseif(length(varargin)==1)
        T=varargin{1};
        n=1000;
    else
        T=varargin{1};
        n=varargin{2};
    end
    z=zeros(n,3);
    e=self.e;
    e=e(2:5);
    flag=true(n,1);
    TRun=zeros(n,1);
    dtotal=1;
    v=[1 0 0
        -1 0 0
        0 1 0
        0 -1 0
        0 0 1
        0 0 -1];
    while any(flag)
        time=TRun(flag);
        num=z(flag,:);
        l=[(dtotal-num(:,1))*t1,num(:,1)*t0,k1*num(:,1),dm*num(:,2),kp*num(:,2),dp*num(:,3)];
        t=exprnd(l);
        [Tspan,id]=min(t,[],2);
        vc=v(id,:);
        time=time+Tspan;
        num=num+vc;
        TRun(flag)=time;
        z(flag,:)=num;
        flag(TRun>=T)=false;
        id=repmat(1:4,n,1);
        zd=repmat(z(:,3),1,4);
        m=mean(zd.^id,1);
        f=abs(m-e)./e;
        if(all(f<0.05))
            break;
        end
    end
    self.simudata=z;
end
