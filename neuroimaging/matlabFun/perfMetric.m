function res = perfMetric(y,yHat,k)
warning('off','stats:perfcurve:SubSampleWithMissingClasses')
averageWR = 1;
if ~exist('y','var')
    res = struct(...
        'y',[],...
        'yHat',[],...
        'nObs',[],...
        'hit',[],...
        'acc',[],...
        'acc_CI5',[],...
        'acc_CI95',[],...
        'acc_thresh',[],...
        'acc_p',[],...);
        'auc',[],...
        'auc_CI5',[],...
        'auc_CI95',[],...
        'distT',[],...
        'distT_p',[],...
        'distMean',[],...
        'distStd',[]);
    return
end

if averageWR
    nRun = length(unique(k))*length(unique(y));
    y = mean(reshape(y,[length(y)/nRun nRun]),1)';
    if nRun~=size(yHat,1)
        error('double-check that')
    end
    yHat = mean(reshape(yHat,[size(yHat,1)/nRun nRun size(yHat,3) size(yHat,4)]),1);
    yHat = permute(yHat,[2 3 4 1]); % run x t x wSect
end

res.y = {y};
res.nObs = length(y);
res.yHat = {yHat};

% Contitnue only if we have a pair of conditions
if length(unique(y))>2
    resTmp = perfMetric;
    fieldList = fields(res);
    for fieldInd = 1:length(fieldList)
        resTmp.(fieldList{fieldInd}) = res.(fieldList{fieldInd});
    end
    res = resTmp;
    return
end

% acc
if size(yHat,2)==1
    res.hit = sum((yHat<0)+1==y);
    res.acc = res.hit./res.nObs;
    [~,pci] = binofit(res.hit,res.nObs,0.1);
    res.acc_CI5 = pci(1);
    res.acc_CI95 = pci(2);
    [~,pci] = binofit(res.nObs/2,res.nObs,0.1);
    res.acc_thresh = pci(2);
    res.acc_p = binocdf(res.hit,res.nObs,0.5,'upper');
else
    res.hit = [];
    res.acc = [];
    res.acc_CI5 = [];
    res.acc_CI95 = [];
    res.acc_thresh = [];
    res.acc_p = [];
end
for tInd = 1:size(yHat,2)
    % auc
    [~,~,~,auc] = perfcurve(y,yHat(:,tInd),1,'NBOOT',2^10);
    res.auc(:,tInd) = auc(1);
    res.auc_CI5(:,tInd) = auc(2);
    res.auc_CI95(:,tInd) = auc(3);
    % distT
    [~,P,~,STATS] = ttest(yHat(y==1,tInd),yHat(y==2,tInd));
    res.distT(:,tInd) = STATS.tstat;
    res.distT_p(:,tInd) = P;
end
% dist
dist = yHat(y==1,:) - yHat(y==2,:);
res.distMean = mean(dist,1);
res.distStd = std(dist,[],1);

