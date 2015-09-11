function self=plots(self,varargin)
%plots is a method of the class FPDAnalyzer which plots the approximated
%probability mass function of the protein level distribution.
%  z=obj.plots() plots the approximated protein distribution defined in
%  obj. All of the plot are automatically set. 
%  z=obj.plots(sim) also plot the histogram of the protein levels in the
%  obj.simudata. If obj.simudata is not empty, the data will be used
%  directly, else obj.Simulater(1000,1000) will be called and the result is
%  used. 
%  z=obj.plots(...,T,n) call obj.Simulater(T,n) to update obj.simudata no
%  matter the histogram is shown.
%  The output is the updated instance. 
    if(isempty(varargin))
        sim=false;
    elseif(length(varargin)==1)
        sim=varargin{1};
        T=1000;
        n=1000;
    elseif(length(varargin)==3)
        sim=varargin{1};
        T=varargin{2};
        n=varargin{3};
    end
    if(isempty(self.plotdata))
        mu=self.e(2);
        sigma=sqrt(self.e(3)-self.e(2)^2);
        if(mu-5*sigma>0)
            x=ceil(mu-5*sigma):ceil(mu+5*sigma);
        else
            x=0:ceil(mu+5*sigma);
        end
        y=self.Analyzer(x);
        self.plotdata=[x(:) y(:)];
    end
    plot(self.plotdata(:,1),self.plotdata(:,2));
    xlabel('Number of protein');
    ylabel('Frequency');
    if(sim)
        if(isempty(self.simudata))
            self=self.Simulater(T,n);
        end
        hold on
        histogram(self.simudata(:,3),'Normalization','pdf');
        legend('EFPD','Gillespie');
    end
end
