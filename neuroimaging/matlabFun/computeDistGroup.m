function dist = computeDistGroup(r,absFalg,tFlag)
if ~exist('absFalg','var')
    absFalg = 0;
end
if ~exist('tFlag','var')
    tFlag = 0;
end

%Extract within-subFold average for each label
yhat1 = nan(sum(r.crossVal(:,1,1)==1&r.y(:,1,1)==1),size(r.yhat,2),size(r.yhat,3),size(r.yhat,4),size(r.yhat,5));
yhat2 = nan(sum(r.crossVal(:,1,1)==1&r.y(:,1,1)==1),size(r.yhat,2),size(r.yhat,3),size(r.yhat,4),size(r.yhat,5));
for subj = 1:size(r.crossVal,3)
    for rep = 1:size(r.crossVal,2)
        crossVal = unique(r.crossVal(:,rep,subj));
        crossVal = crossVal(~isnan(crossVal));
        for crossValInd = 1:length(crossVal)
            curCrossVal = r.crossVal(:,rep,subj)==crossVal(crossValInd);
            curLabel = r.y(:,rep,subj)==1;
            yhat1(:,:,crossVal(crossValInd),rep,subj) = (r.yhat(curCrossVal&curLabel,:,crossVal(crossValInd),rep,subj));
            curLabel = r.y(:,rep,subj)==2;
            yhat2(:,:,crossVal(crossValInd),rep,subj) = (r.yhat(curCrossVal&curLabel,:,crossVal(crossValInd),rep,subj));
        end
    end
end
%Compute distance
if tFlag
     [~,~,~,STATS] = ttest2(yhat1,yhat2);
     dist = shiftdim(STATS.tstat,1);
else
    dist = shiftdim(mean(yhat1 - yhat2,1),1);
end

%Absolute value to avoid anti-learning effects
if absFalg
    dist = abs(dist);
end
%averatge across subfolds
dist = shiftdim(nanmean(dist,1),1);
%averatge across folds
dist = shiftdim(nanmean(dist,1),1);


