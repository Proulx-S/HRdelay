function dist = computeDist(r,absFalg,tFlag)
if ~exist('absFalg','var')
    absFalg = 0;
end
if ~exist('tFlag','var')
    tFlag = 0;
end

%Extract within-subFold average for each label
for rep = 1:size(r.crossVal,2)
    crossVal = unique(r.crossVal(:,rep));
    for crossValInd = 1:length(crossVal)
        curCrossVal = r.crossVal(:,rep)==crossVal(crossValInd);
        curLabel = r.y(:,rep)==1;
        yhat1(:,:,crossVal(crossValInd),rep) = (r.yhat(curCrossVal&curLabel,:,crossVal(crossValInd),rep));
        curLabel = r.y(:,rep)==2;
        yhat2(:,:,crossVal(crossValInd),rep) = (r.yhat(curCrossVal&curLabel,:,crossVal(crossValInd),rep));
    end
end
%Compute distance
if tFlag
     [~,~,~,STATS] = ttest2(yhat1,yhat2);
     dist = squeeze(STATS.tstat);
else
    dist = squeeze(mean(yhat1 - yhat2,1));
end

%Absolute value to avoid anti-learning effects
if absFalg
    dist = abs(dist);
end
%averatge across subfolds
dist = squeeze(mean(dist,1));
%averatge across folds
dist = squeeze(mean(dist,1));


