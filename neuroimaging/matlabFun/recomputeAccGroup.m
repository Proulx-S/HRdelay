function acc = recomputeAccGroup(r,rectifyAcc)
if ~exist('rectifyAcc','var')
    rectifyAcc = 0;
end
if ~rectifyAcc
    %collapse folds
    yhat = nanmean(r.yhat,3);
    %averatge across subfolds
    yhat = squeeze(nanmean(yhat,2));
    %recompute acc
    for subj = 1:size(r.y,3)
        for rep = 1:size(r.y,2)
            acc(1,rep,subj) = (sum(yhat(r.y(:,rep,subj)==1,rep,subj)>=0,1)+sum(yhat(r.y(:,rep,subj)==2,rep,subj)<0,1))./size(r.y(~isnan(r.y(:,rep,subj)),rep,subj),1)*100;
        end
    end
else
    %recompute acc
    for subj = 1:size(r.yhat,5)
        for rep = 1:size(r.yhat,4)
            for fold = 1:size(r.yhat,3)
                for subfold = 1:size(r.yhat,2)
                    curInd = ~isnan(r.yhat(:,subfold,fold,rep,subj));
                    acc(subfold,fold,rep,subj) = (sum(r.yhat(curInd&r.y(:,rep,subj)==1,subfold,fold,rep,subj)>=0) + sum(r.yhat(curInd&r.y(:,rep,subj)==2,subfold,fold,rep,subj)<0)) / sum(curInd);
                    acc(subfold,fold,rep,subj) = abs(acc(subfold,fold,rep,subj)-0.5)+0.5;
                end
            end
        end
    end
    acc = squeeze(nanmean(nanmean(acc,1),2)) * 100;
end