function r = averageFoldsTe(r)
if isfield(r,'acc')
    r.acc = mean(r.acc,2);
end
if isfield(r,'w')
    r.w= mean(r.w,1);
    r.a= mean(r.a,1);
end
for i = 1:size(r.yhat,2)
    r.yhat(r.crossVal==i,i)
    for ii = 1:length(unique(r.crossVal))
    r.yhat(r.crossVal==ii,1)
    end
r.yhat = nanmean(r.yhat,2);
if isfield(r,'pred')
    r.pred = nanmean(r.pred,2);
end

if isfield(r,'subLevel')
    for subInd = 1:length(r.subLevel)
        r.subLevel{subInd} = averageFolds(r.subLevel{subInd});
    end
end


