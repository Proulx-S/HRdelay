function metric = applyPhaseLockingThresh(metric)

phaseLockInd = zeros(size(metric.delayRect,1),2);
for ii = 1:2
    tmp     = sec2rad(metric.delayRect(:,:,ii),metric.fittingParam.period);
    pval    = nan(size(tmp,1),1);
    z       = nan(size(tmp,1),1);
    for i = 1:size(tmp,1)
        [pval(i), z(i)] = circ_rtest(tmp(i,:));
    end
    phaseLockInd(:,ii) = pval<0.05;
end
% length(find(diff(phaseLockInd,[],2)))/numel(phaseLockInd)
% length(find(diff(phaseLockInd,[],2)))/size(phaseLockInd,1)
% length(find(any(phaseLockInd,2)))/size(phaseLockInd,1)
% length(find(all(phaseLockInd,2)))/length(find(any(phaseLockInd,2)))

phaseLockInd = any(phaseLockInd,2);


for i = 1:length(metric.dataFields)
    curField = metric.dataFields{i};
    metric.(curField) = metric.(curField)(phaseLockInd,:,:);    
end
metric.voxInd = metric.voxInd(phaseLockInd,:);