classdef EFPD
    properties
        para;
        order;
        e;
        simudata=[];
    end
    
    properties (SetAccess=private,GetAccess=private)
        ecount;
        plotdata=[];
        fun;
    end
    
    methods
        function self=EFPD(para,order)
%EFPD is a MATLAB class who realizes the EFPD algorithm in 
%approximating the protein level distribution at steady state of the two 
%state gene expression model when the parameters in the model are offered. 
%
%PROPERTIES:
%'para': the array of parameters in the two state geneexpression model. 
%para=[t_0,t_1,k_m,k_p,d_m,d_p]. All six elements must be nonnegative real 
%numbers and must be given when an instance of the class is constructed. 
%'order': the order of moments used in EFPD approximation to the protein 
%level distribution. It should be a integer with 4<=order<=10 and must be 
%given when an instance of the class is constructed. 
%"e": the moments of the protein distribution at steady state calculated by
%EFPD. length(e)=order+1 and e(1)=1. It will be calculated automatically 
%when an instance of the class is constructed. It should not be assigned 
%or modified manually.
%'simudata': the simulated data constructed by the Gillespie algorithm. 
%When an instance is constructed it is [] and could be constructed by the 
%Simulater method. It is an n-by-3 matrix, in which, the first column is
%the state of the gene (1: active, 0: inactive), the second column the
%number of mRNA and the third column the number of protein. 
%
%METHODS:
%'Analyzer':
%Analyzer is a method of the class EFPD which calculates the
%approximated probability of given protein level(s).
%  z=obj.Analyzer(x) returns the approximated probability of protein
%  level(s) given in the one dimensional nonnegative integer array x for
%  the model defined in obj. 
%'plots':
%plots is a method of the class EFPD which plots the approximated
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
%'Simulater':
%Simulater is a method of the class EFPD which using the Gillespie
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
%
%USAGE:
%Use obj=EFPD(para,order) to initialize an instance of the class.
%para=[t_0,t_1,k_m,k_p,d_m,d_p] where all the parameters are nonegative
%real numbers. order should be an nonnegative integer between 4 and 10.

            self.para=para;
            self.order=order;
            self.e=self.varnewfindmoments();
            self.ecount=self.e;
            for n=1:order
                s=0;
                for k=1:(n+1)
                    s=s+nchoosek(n+1,k)*self.e(n+2-k);
                end
                self.ecount(n+1)=s/(n+1);
            end
            mu=self.ecount(2);
            sigma=sqrt(self.ecount(3)-self.ecount(2)^2);
            if(mu-5*sigma>0)
                self.fun=@(x)self.NewHermiteApproximate(x);
            else
                self.fun=@(x)self.LaguerreApproximate(x);
            end
        end
        z=Analyzer(self,x);
        self=Simulater(self,varargin);
        self=plots(self,varargin);
    end
    
    methods (Access=private)
        z=varnewfindmoments(self);
        z=NewHermiteApproximate(self,x);
        z=LaguerreApproximate(self,x);
    end
end
