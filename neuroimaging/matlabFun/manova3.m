function [stat,PVAL] = manova3(Xmat,Ymat,C,A,D,statLabel)
[Beta,DFE,Cov] = fitrm2(Xmat,Ymat);
SSE = DFE * Cov;

stat = nan(size(A,1),1);
PVAL = nan(size(A,1),1);
i = 0;
for withinTestInd = 1:length(C)
    for betweenTestInd = 1:length(A)
        i = i+1;
        [stat(i),PVAL(i),~] = getStats(Xmat,A{betweenTestInd},Beta,C{withinTestInd},D,SSE,statLabel);
    end
end