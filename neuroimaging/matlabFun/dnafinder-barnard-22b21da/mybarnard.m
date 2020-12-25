function STATS=mybarnard(x,varargin)
%BARNARD Barnard's Exact Probability Test.
%There are two fundamentally different exact tests for comparing the equality
%of two binomial probabilities – Fisher’s exact test (Fisher, 1925), and
%Barnard’s exact test (Barnard, 1945). Fisher’s exact test (Fisher, 1925) is
%the more popular of the two. In fact, Fisher was bitterly critical of
%Barnard’s proposal for esoteric reasons that we will not go into here. 
%For 2 × 2 tables, Barnard’s test is more powerful than Fisher’s, as Barnard
%noted in his 1945 paper, much to Fisher’s chagrin. Anyway, perhaps due to its
%computational difficulty the Barnard's is not widely used. This function is
%completely vectorized and without for...end loops, and so, the computation is
%very fast.
%The Barnard's exact test is a unconditioned test for it generates the exact
%distribution of the Wald statistic T(X),
%
%           T(X) = abs((p(a) - p(b))/sqrt(p*(1-p)*((1/c1)+(1/c2)))),
%   where,
%           p(a) = a/c1, p(b) = b/c2 and p = (a+b)/n, 
%
%by considering all tables X and calculates P(np) for all possible values of 
%np€(0,1). 
%Under H0, the probability of observing any generic table X is
%
%            /Cs1\  /Cs2\   (I+J)       [N-(I+J)]
% P(X|np) =  |    | |    |*np     *(1-np)
%            \ I /  \ J /
%
%Then, for any given np, the exact p-value of the observed Table Xo is 
%        __
%        \
%  p(np)=/_P(X|np)
%        T(X)>=T(Xo)
%
%Barnard suggested that we calculate p(np) for all possible values of np€(0,1)
%and choose the value, np*, say, that maximizes p(np): PB=sup{p(np): np€(0,1)}.
%
%   Syntax: function [STATS]=mybarnard(x,plts,Tbx) 
%      
%      
%     Inputs:
%           X - 2x2 data matrix
%           PLTS - Flag to set if you don't want (0) or want (1) view the plots (default=0)
%           Tbx - is the granularity of the np array (how many points in the
%           interval (0,1) must be considered to determine np* (default=100).
%     Output:
%         A table with:
%         - Wald statistic, Nuisance parameter and P-value
%         - Plot of the nuisance parameter PI against the corresponding P-value for
%           all the PI in (0, 1). It shows the maximized PI where it attains the
%           P-value.
%        If STATS nargout was specified the results will be stored in the STATS
%        struct.
%
%   Example:
%
%                                    Vaccine
%                               Yes           No
%                            ---------------------
%                    Yes         7            12
% Infectious status                 
%                     No         8            3
%                            ---------------------
%                                       
%   Calling on Matlab the function: 
%             mybarnard([7 12; 8 3])
%
%   Answer is:
%
%     Tables      Size      Wald_stat    Nuisance    one_tailed_p_value    two_tailed_p_value
%     ______    ________    _________    ________    __________________    __________________
% 
%     100       16    16    1.8943       0.66663     0.034074              0.068148   
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2009) MyBarnard: a very compact routine for Barnard's exact test on 2x2 matrix
% http://www.mathworks.com/matlabcentral/fileexchange/25760

%Input Error handling
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'real','finite','integer','nonnegative','nonnan','size',[2 2]}));
addOptional(p,'plts',0, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
addOptional(p,'Tbx',100, @(x) validateattributes(x,{'numeric'},{'real','finite','integer','positive','nonnan','scalar'}));
parse(p,x,varargin{:});
Tbx=p.Results.Tbx; plts=p.Results.plts;
clear p

Cs=sum(x); N=sum(Cs); %Columns sum and total observed elements
%Rearrange the matrix if necessary.
flip=0;
if ~issorted(Cs)
     x=fliplr(x);
     Cs=sort(Cs);
     flip=1;
end

%First of all, we must compute the Wald's statistic.
%           T(X) = (p(a) - p(b))/sqrt(p*(1-p)*((1/Cs1)+(1/Cs2))),
%   where:
% I can vary from 0 to Cs(1)
% J can vary from 0 to Cs(2)
% p(a) = I/Cs1, p(b) = J/Cs2 and p = (I+J)/N
% It is possible to compute all the T(X) values, without using the for...end loops

[I,J]=ndgrid(0:1:Cs(1), 0:1:Cs(2));

% I and J will be 2 square matrix
%I=  0    0    0  ... ...  0
%    1    1    1  ... ...  1
%   ...  ...  ... ... ... ...
%   ...  ...  ... ... ... ...
%   Cs1  Cs1  Cs1 ... ... Cs1
%
%J=  0    1    2  ... ...  Cs2
%    0    1    2  ... ...  Cs2
%   ...  ...  ... ... ...  ...
%   ...  ...  ... ... ...  ...
%    0    1    2  ... ...  Cs2
%
%TX will be a Cs(1)+1 x Cs(2)+1 matrix.
%For two I and J combinations, T(X) shows a singularity:
%i.e. I=0 J=0 T(X)=0/0
%     I=Cs(1) J=Cs(2) T(X)=0/0
%So, for these two points T(X) will be set to 0

TX=(I./Cs(1)-J./Cs(2))./realsqrt(((I+J)./N).*(1-((I+J)./N)).*sum(1./Cs)); %computes the Wald's statistics
TX([1 end])=0; %resolve the singularities
TX0=abs(TX(x(1)+1,x(3)+1)); %catch the observed Table (TXo), taking the 0 in account
idx=TX>=TX0; %set the index of all T(X)>=TXo

%Idx will be a Cs(1)+1 x Cs(2)+1 matrix where cells are 1 if TX>=TXo otherwise 0

%Now we must introduce the nuisance parameter.
%Under H0, the probability of observing any generic table X is
%
%            /Cs1\  /Cs2\   (I+J)       [N-(I+J)]
% P(X|np) =  |    | |    |*np     *(1-np)
%            \ I /  \ J /
%
%Then, for any given np, the exact p-value of the observed Table Xo is 
%        __
%        \
%  p(np)=/_P(X|np)
%        T(X)>=T(Xo)
%What shall we do with the unknown nuisance parameter np in the above p-value?
%Barnard suggested that we calculate p(np) for all possible values of np€(0,1)
%and choose the value, np*, say, that maximizes p(np).
% PB=sup{p(np): np€(0,1)}
%
%Expanding P(X|np) we have:
%
%                  (I+J)        [N-(I+J)]
%  Cs1! * Cs2! * np      * (1-np)
%  --------------------------------------
%     I! * (Cs1-I)! * J! * (Cs2-J)!
%
%The factorials growth very quickly, so to avoid and overflow:
%x!=prod(1:x) => log(x!)=sum(log(1:x))
%x!=G(x+1) => log(x!)=Glog(x+1) where G is the Gamma function and Glog is the
%logarithm of the Gamma function.
%Using the logarithm, multiplications and divisions became sums (and computers
%love sums much more than multiplications...)
%log(a*b)=log(a)+log(b); log(a/b)=log(a)-log(b); log(a^b)=b*log(a)
%So:
%log(P(X|np))=[Glog(Cs1+1)+Glog(Cs2+1)-Glog(I+1)-Glog(Cs1-I+1)-Glog(J+1)-Glog(Cs2-J+1)]+[(I+J)*log(np)]+{[N-(I+J)]*log(1-np)]}
% 
%Thus, to compute log(P(X|np)) we should implement three nested for...end loops
%for np=0:step:1
%  for I=0:Cs1
%      for J=0:Cs2
%          compute S=log(P(np,I,J))
%      end
%  end
%  P(np)=sum(exp(S(TX>=TXo)));
% end
%p-value=Max(P); nuisance parameter=np(P==p-value)
%
%This approach is very time consuming and is not efficient. Consider that is 
%the users that establishes the duration of the first loop: little step values
%(Tbx variable) raise up the accuracy of the nuisance parameter computations but
%the time of computing becames highest. Moreover, in Matlab you can do the same
%without using the for...end loops. 
%First of all, the quantity between the first square brackets depends only on I
%and J but not on np
%log(P(X|np))=CF+[(I+J)*log(np)]+{[N-(I+J)]*log(1-np)]}
%I and J were implemented using ndgrid function and so CF, as TX, will be a
%(Cs1+1) x (Cs2+2) matrix. But, also [(I+J)*log(np)] and {[N-(I+J)]*log(1-np)]}
%will be (Cs1+1) x (Cs2+2) matrix. These matrix will be summed Tbx times...

%Set the basic parameters...
A=[1 1 Tbx]; B=Cs+1; npa = linspace(0.0001,0.9999,Tbx); 
LP=log(npa); ALP=log(1-npa); E=repmat(I+J,A); F=N-E;

clear Cs N 

%Generate a 3D matrix (Cs1+1 x Cs2+1 x Tbx): this is a "box" where each of Tbx 2D
%slice is the Cs1+1 x Cs2+1 matrix CF
CF=repmat(sum(gammaln(B))-(gammaln(I+1)+gammaln(J+1)+gammaln(B(1)-I)+gammaln(B(2)-J)),A);

clear I J 
%Now we have two 1xTbx arrays (LP and ALP) and two Cs1+1 x Cs2+1 x Tbx "boxes"
%(E and F). As Peter J. Acklam wrote in his "Matlab array manipulation tips and
%tricks" (http://home.online.no/~pjacklam/matlab/doc/mtt/doc/mtt.pdf), there is
%a no for...loop solution to multiply each 2D slice with the corresponding
%element of a vector. 
%Assume X is an m-by-n-by-p array and v is a row vector with length p. 
%How does one write:
%
%  Y = zeros(m, n, p);
%  for i = 1:p
%     Y(:,:,i) = X(:,:,i) * v(i);
%  end
%
%with no for-loop? One way is to use:
%  Y = X .* repmat(reshape(v, [1 1 p]), [m n]);
% 
%So we can construct the two others 3D boxes:
%Box1=CF
%Box2=E.*repmat(reshape(LP,[1 1 Tbx]),[Cs1+1 Cs2+1])=E.*repmat(reshape(LP,A),B)
%Box3=F.*repmat(reshape(ALP,[1 1 Tbx]),[Cs1+1 Cs2+1])=F.*repmat(reshape(ALP,A),B)
%Finally, we can obtain the S box that is the sum of these three boxes (remember
%that we are using logarithms, so we must convert them using the exp function)
S=exp(CF+E.*repmat(reshape(LP,A),B)+F.*repmat(reshape(ALP,A),B));

clear CF LP ALP F E
%Now we are at the last point: for each 2D slice of the S box we must sum the
%values corresponding to TX>=TXo, obtaining a new vector P.
%To use the logical indexing tecnique, we must replicate the idx matrix:
%S(repmat(idx,A)) is a column vector that has L*Tbx row; 
%L=sum(idx(idx==1)) is the number of 1 in the idx matrix.
%Using the reshape function we can obtain a new LxTbx matrix: in each column
%we'll have the P(X|np & TX>=Txo) and using the function sum we'll, finally,
%obtain the P vector.
P=sum(reshape(S(repmat(idx, A)),sum(idx(idx==1)),Tbx));

clear A idx 

PV=max(P); %The p-value is tha max value of the P vector;
np=npa(P==PV); %The nuisance parameter is the value of the npa array coinciding with PMax

%display results
disp(table(Tbx,B,TX0,np,PV,min(2*PV,1),...
    'VariableNames',{'Tables','Size','Wald_stat','Nuisance','one_tailed_p_value','two_tailed_p_value'}));
if nargout
    STATS.TX0=TX0;
    STATS.p_value=PV;
    STATS.nuisance=np;
end
clear B TX0

if plts
    %display plot
    figure('Color',[1 1 1],'outerposition',get(groot,'ScreenSize'));
    subplot(1,2,1)
    hold on; 
    patch([0 npa 1],[0 P 0],'b');
    plot([0 np],[PV PV],'k--'); 
    plot([np np],[0 PV],'w--'); 
    plot(np,PV,'ro','MarkerFaceColor','red'); 
    hold off
    title('Barnard''s exact p-value as a function of nuisance parameter np.','FontName','Arial','FontSize',14,'FontWeight','Bold');
    xlabel('Nuisance parameter (np)','FontName','Arial','FontSize',14,'FontWeight','Bold');
    ylabel('p-value  [p(TX >= TXO | np)]','FontName','Arial','FontSize',14,'FontWeight','Bold');
    axis square
    
    
    txt=['P[T(X)|pi=' num2str(np) ']'];
    clear npa np
    A=TX(:);
    B=reshape(S(:,:,P==PV),numel(TX),1);
    clear S P PV TX
    [As,idx]=sort(A);
    Bs=B(idx);
    clear A B idx S P PV TX
    Xplot=unique(As);
    Yplot=zeros(size(Xplot));
    for I=1:length(Xplot)
        Yplot(I)=sum(Bs(As==Xplot(I)));
    end
    clear As Bs I
    subplot(1,2,2)
    bar(Yplot) %bar(Xplot,Yplot)
    axis square tight
    xt=get(gca,'Xtick'); xtl=cell(1,length(xt));
    for I=1:length(xt)
        xtl{I}=sprintf('%0.1f',Xplot(xt(I)));
    end
    set(gca,'Xtick',xt,'XTickLabel',xtl)
    clear I xt*
    title('Distribution of Wald Statistic for Barnard''s Exact test','FontName','Arial','FontSize',12,'FontWeight','Bold');
    xlabel('T(X)','FontName','Arial','FontSize',12,'FontWeight','Bold');
    ylabel(txt,'FontSize',12,'FontWeight','Bold');
end