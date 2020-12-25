# barnard
Barnard's Exact Probability Test.<br/>
There are two fundamentally different exact tests for comparing the equality
of two binomial probabilities – Fisher’s exact test (Fisher, 1925), and
Barnard’s exact test (Barnard, 1945). Fisher’s exact test (Fisher, 1925) is
the more popular of the two. In fact, Fisher was bitterly critical of
Barnard’s proposal for esoteric reasons that we will not go into here. 
For 2 × 2 tables, Barnard’s test is more powerful than Fisher’s, as Barnard
noted in his 1945 paper, much to Fisher’s chagrin. Anyway, perhaps due to its
computational difficulty the Barnard's is not widely used. This function is
completely vectorized and without for...end loops, and so, the computation is
very fast.
The Barnard's exact test is a unconditioned test for it generates the exact
distribution of the Wald statistic T(X),

          T(X) = abs((p(a) - p(b))/sqrt(p*(1-p)*((1/c1)+(1/c2)))),
  where,
          p(a) = a/c1, p(b) = b/c2 and p = (a+b)/n, 

by considering all tables X and calculates P(np) for all possible values of 
np€(0,1). 
Under H0, the probability of observing any generic table X is

           /Cs1\  /Cs2\   (I+J)       [N-(I+J)]
P(X|np) =  |    | |    |*np     *(1-np)
           \ I /  \ J /

Then, for any given np, the exact p-value of the observed Table Xo is 
       __
       \
 p(np)=/_P(X|np)
       T(X)>=T(Xo)

Barnard suggested that we calculate p(np) for all possible values of np€(0,1)
and choose the value, np*, say, that maximizes p(np): PB=sup{p(np): np€(0,1)}.

  Syntax: function [STATS]=mybarnard(x,plts,Tbx) 
     
     
    Inputs:
          X - 2x2 data matrix
          PLTS - Flag to set if you don't want (0) or want (1) view the plots (default=0)
          Tbx - is the granularity of the np array (how many points in the
          interval (0,1) must be considered to determine np* (default=100).
    Output:
        A table with:
        - Wald statistic, Nuisance parameter and P-value
        - Plot of the nuisance parameter PI against the corresponding P-value for
          all the PI in (0, 1). It shows the maximized PI where it attains the
          P-value.
       If STATS nargout was specified the results will be stored in the STATS
       struct.

  Example:

                                   Vaccine
                              Yes           No
                           ---------------------
                   Yes         7            12
Infectious status                 
                    No         8            3
                           ---------------------
                                      
  Calling on Matlab the function: 
            mybarnard([7 12; 8 3])

  Answer is:

    Tables      Size      Wald_stat    Nuisance    one_tailed_p_value    two_tailed_p_value
    ______    ________    _________    ________    __________________    __________________

    100       16    16    1.8943       0.66663     0.034074              0.068148   

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2009) MyBarnard: a very compact routine for Barnard's exact test on 2x2 matrix
http://www.mathworks.com/matlabcentral/fileexchange/25760
