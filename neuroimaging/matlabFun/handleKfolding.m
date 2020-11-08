function [curd,tmpp] = handleKfolding(d,p)
curd = d;
tmpp = p;

if ischar(p.k) && p.trialsPerObs==8
    p.k = p.nObs/2;
end
if mod(length(d.label)/p.k,1)
    error('number of data point is not a multiple of k (programmer in a rush)')
end

curd.crossVal = getRandCrossValFolds(d,p);
