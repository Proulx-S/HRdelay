function acc = recomputeAcc(r,rectifyAcc)
if ~exist('rectifyAcc','var')
    rectifyAcc = 0;
end
if ~rectifyAcc
    %collapse folds
    yhat = nanmean(r.yhat,3);
    %averatge across subfolds
    yhat = squeeze(mean(yhat,2));
    %recompute acc
    for rep = 1:size(r.y,2)
        acc(1,rep) = (sum(yhat(r.y(:,rep)==1,rep)>=0,1)+sum(yhat(r.y(:,rep)==2,rep)<0,1))./size(r.y,1)*100;
    end
else
    %recompute acc
    for rep = 1:size(r.yhat,4)
        for fold = 1:size(r.yhat,3)
            for subfold = 1:size(r.yhat,2)
                curInd = ~isnan(r.yhat(:,subfold,fold,rep));
                acc(subfold,fold,rep) = (sum(r.yhat(curInd&r.y(:,rep)==1,subfold,fold,rep)>=0) + sum(r.yhat(curInd&r.y(:,rep)==2,subfold,fold,rep)<0)) / sum(curInd);
                acc(subfold,fold,rep) = abs(acc(subfold,fold,rep)-0.5)+0.5;
            end
        end
    end
    acc = squeeze(mean(mean(acc,1),2))' * 100;
end