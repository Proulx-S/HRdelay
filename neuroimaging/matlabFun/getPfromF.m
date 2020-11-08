function [FDR,pGood] = getPfromF(F,df)
%Convert F to p
%Convert to p and Z, with trick to avoid p=0 and z=+inf at very large F
p = fcdf(F,df.pFull-df.pReduced,df.n-df.pFull);
q = fcdf(1./F,df.n-df.pFull,df.pFull-df.pReduced);
pGood = nan(size(p));
pGood(p<=0.5) = 1-p(p<=0.5);
pGood(p>0.5) = q(p>0.5);
%         [FDR,Q,PIO] = mafdr(pGood);
FDR = mafdr(pGood(:),'BHFDR',true);
tmp = nan(size(pGood));
tmp(1:numel(pGood)) = FDR;
FDR = tmp; clear tmp
