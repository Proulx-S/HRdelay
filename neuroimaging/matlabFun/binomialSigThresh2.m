function hitRate = binomialSigThresh2(svm,alpha)
if ~exist('alpha','var')
    alpha = 0.05;
end

hitRate = squeeze(nanmean(svm.r.hitRate,1));
negThresh = binoinv(alpha,svm.p.nObs,0.5)/svm.p.nObs;
posThresh = binoinv(1-alpha,svm.p.nObs,0.5)/svm.p.nObs;
hitRate(hitRate>negThresh & hitRate<posThresh) = nan;
