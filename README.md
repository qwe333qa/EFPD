# EFPD
----------------------------------------------------------
README FILE OF THE EFPD APPLICATION IN MATLAB
----------------------------------------------------------

DESCRIPTION

This program package is the implementation of the EFPD algorithm in MATLAB. This program is packed into a MATLAB class. EFPD is an algorithm which constructs an approximation of the protein level distribution of the two state gene expression model. The model could be shown as:

         k_m        k_p
 Active------>mRNA------->Protein
 |    /\       |            |
 |t_0 |t_1     |d_m         |d_p
 |    |        |            |
 \/   |        \/           \/
Inactive       0            0

The model is charactered by the six-tuple [t_0,t_1,k_m,k_p,d_m,d_p], in which, the gene is switched between inactive and active state with rates t_0 (active to inactive) and t_1 (inactive to active) respectively. when gene is active, mRNAs are transcribed with a rate k_m. Proteined are translated with a rate k_p. mRNAs and proteins are degradated with rates d_m and d_p respectively. All the reactions are independent with each other.  

The EFPD takes the six-tuple of the parameters as the input, then calculates the lower order moments of the protein distribution. Due to the linearity of the model, the calculation of the moments is straightforward. Fourier theories are used to construct an approximation of the distribution. 

FILE STRUCTURE

@EFPD  %The folder cantaining the class.
|--Analyzer.m
|--EFPD  The class definition.
|--LaguerreApproximate.m
|--NewHermiteApproximate.m
|--plots.m
|--README.txt
|--Simulater.m
|--varnewfindmoments.m

PROPERTIES

para: the array of parameters in the two state geneexpression model. para=[t_0,t_1,k_m,k_p,d_m,d_p]. All six elements must be nonnegative real numbers and must be given when an instance of the class is constructed. 

order: the order of moments used in EFPD approximation to the protein level distribution. It should be a integer with 4<=order<=10 and must be given when an instance of the class is constructed. 

e: the moments of the protein distribution at steady state calculated by EFPD. length(e)=order+1 and e(1)=1. It will be calculated automatically when an instance of the class is constructed. It should not be assigned or modified manually.

simudata: the simulated data constructed by the Gillespie algorithm. When an instance is constructed it is [] and could be constructed by the Simulater method. It is an n-by-3 matrix, in which, the first column is the state of the gene (1: active, 0: inactive), the second column the number of mRNA and the third column the number of protein. 

PRIVATE PROPERTIES
These private properties can't be reached outside the class.

ecount: moments of the continuous version of the distribution approximated. It's calculated by e.

plotdata: the data used in plots. It is assigned by the plots method. 

fun: the function handle in which the approximator used is packed. 

METHODS

Analyzer:
Analyzer is a method of the class EFPD which calculates the approximated probability of given protein level(s).
  z=obj.Analyzer(x) returns the approximated probability of protein
  level(s) given in the one dimensional nonnegative integer array x for
  the model defined in obj. 

plots:
plots is a method of the class EFPD which plots the approximated probability mass function of the protein level distribution.
  z=obj.plots() plots the approximated protein distribution defined in
  obj. All of the plot are automatically set. 
  z=obj.plots(sim) also plot the histogram of the protein levels in the
  obj.simudata. If obj.simudata is not empty, the data will be used
  directly, else obj.Simulater(1000,1000) will be called and the result is
  used. 
  z=obj.plots(...,T,n) call obj.Simulater(T,n) to update obj.simudata no
  matter the histogram is shown.
  The output is the updated instance. 

Simulater:
Simulater is a method of the class EFPD which using the Gillespie algorithm to simulate the reaction path of gene expression from a total inactive state (gene inactive, no mRNA and protein).
  z=obj.Simulater(T,n) simulate a reaction patt of n of the model defined
  in obj. The program will stop when time T is reached or all the first
  four moments of the simulated protein population are differ no more than
  5 from the theoretical ones. the property simudata is updated by the
  newly generated data.
  z=obj.Simulater(T) take n=1000 as the default value. 
  z=obj.Simulater() also take T=1000 as default value.
  The output stores the updated instance. 
  The simulation may take a lot of time.

PRIVATE METHODS
These private methods can't be called outside the class.

varnewfindmoments:
The moments calculater. Used in the calculation of moments of protein distribution up to the order given in order with parameters given in para.

NewHermiteApproximate:
The approximater with a normal kernel. 

LaguerreApproximate:
The approximater with a Gamma kernel. 

INSTALL AND USAGE

The package is totally open source and installation-free. To use it, you only need to copy all the files to a directory named @EFPD whose parent directory is on the MATLAB path.

To begin, you need to create an object of the class using the command 
  x=EFPD(para,order);
in matlab. Note the order of parameters passed to "para" and the range of the vlue passed to "order". Then call the methods or properties of x on your need.

Please cite "An efficient and assumption free method to approximate protein level distribution in the two states gene expression model" if you use the package in your research work.
