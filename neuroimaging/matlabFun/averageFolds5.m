function r = averageFolds5(r)
if isfield(r,'acc')
    r.acc = mean(r.acc,2);
end
if isfield(r,'w')
    r.w= mean(r.w,1);
end
if isfield(r,'a')
    r.a= mean(r.a,1);
end
r.yhat = nanmean(r.yhat,2);
if isfield(r,'pred')
    r.pred = nanmean(r.pred,2);
end

if isfield(r,'subLevel')
    for subInd = 1:length(r.subLevel)
        r.subLevel{subInd} = averageFolds5(r.subLevel{subInd});
    end
end


