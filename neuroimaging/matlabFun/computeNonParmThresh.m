function threshAcc = computeNonParmThresh(hitRate,threshPerc)
threshAcc = nan(length(threshPerc),size(hitRate,2));
for curFeat = 1:size(threshAcc,2)
    curHitRate = hitRate(:,curFeat);
    curHitRate = curHitRate(~isnan(curHitRate));
    if ~isempty(curHitRate)
        threshAcc(:,curFeat) = prctile(curHitRate,threshPerc);
    end
end
    