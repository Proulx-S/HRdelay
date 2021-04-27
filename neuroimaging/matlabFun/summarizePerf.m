function [resSess,resSubj,resGroup] = summarizePerf(res_sess)
allField = fields(res_sess);
for i = 1:length(allField)
    if isnumeric(res_sess(1).(allField{i}))
        tLength = length(res_sess(1).(allField{i}));
        if tLength==1
            resSess.(allField{i}) = nan(size(res_sess));
        elseif tLength==0
            resSess.(allField{i}) = [];
        else % for the special case where we compute stats at each time delay relative to stimulus onset
            resSess.(allField{i}) = nan([tLength size(res_sess)]);
        end
        if ~isempty(res_sess(1).(allField{i}))
            resSess.(allField{i})(:) = [res_sess.(allField{i})];
        end
        if length(size(resSess.(allField{i})))>2
            resSess.(allField{i}) = permute(resSess.(allField{i}),[2 3 1]);
        end
    elseif iscell(res_sess(1).(allField{i}))
        resSess.(allField{i}) = cell(size(res_sess));
        resSess.(allField{i})(:) = [res_sess.(allField{i})];
    elseif ischar(res_sess(1).(allField{i}))
        resSess.(allField{i}) = cell(size(res_sess));
        resSess.(allField{i})(:) = {res_sess.(allField{i})};
    elseif isstruct(res_sess(1).(allField{i}))
        resSess.(allField{i}) = cell(size(res_sess));
        resSess.(allField{i})(:) = {res_sess.(allField{i})};
    else
        error('code that')
    end
end
resSess.info = 'subj x sess';
resSess.distT_fdr = resSess.distT_p;
resSess.distT_fdr(:) = mafdr(resSess.distT_p(:),'BHFDR',true);
% resSess = orderfields(resSess,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 23 16 17 18 19 20 21 22]);

% Continue only if having a pair of condtions
if isempty(resSess.acc)
    resSubj  = [];
    resGroup = [];
    return
end

resSubj.y = cell(size(resSess.y,1),1);
resSubj.yHat = cell(size(resSess.yHat,1),1);
for subjInd = 1:size(resSess.y,1)
    resSubj.y{subjInd} = cell2mat(resSess.y(subjInd,:)');
    resSubj.yHat{subjInd} = cell2mat(resSess.yHat(subjInd,:)');
end
resSubj.nObs = sum(resSess.nObs,2);
if ~isempty(resSess.hit)
    resSubj.hit = sum(resSess.hit,2);
    resSubj.acc = resSubj.hit./resSubj.nObs;
    [~,pci] = binofit(resSubj.hit,resSubj.nObs,0.1);
    resSubj.acc_CI5 = pci(:,1);
    resSubj.acc_CI95 = pci(:,2);
    [~,pci] = binofit(resSubj.nObs/2,resSubj.nObs,0.1);
    resSubj.acc_thresh = pci(:,2);
    resSubj.acc_p = binocdf(resSubj.hit,resSubj.nObs,0.5,'upper');
    resSubj.acc_fdr = mafdr(resSubj.acc_p,'BHFDR',true);
else
    resSubj.hit = [];
    resSubj.acc = [];
    resSubj.acc_CI5 = [];
    resSubj.acc_CI95 = [];
    resSubj.acc_thresh = [];
    resSubj.acc_p = [];
    resSubj.acc_fdr = [];
end
resSubj.auc = nan([size(resSubj.y) size(resSubj.yHat{1},2)]);
resSubj.distT = nan([size(resSubj.y) size(resSubj.yHat{1},2)]);
resSubj.distT_p = nan([size(resSubj.y) size(resSubj.yHat{1},2)]);
resSubj.distMean = nan([size(resSubj.y) size(resSubj.yHat{1},2)]);
resSubj.distStd = nan([size(resSubj.y) size(resSubj.yHat{1},2)]);
for subjInd = 1:size(resSubj.y,1)
    for tInd = 1:size(resSubj.yHat{subjInd},2)
        [~,~,~,auc] = perfcurve(resSubj.y{subjInd},resSubj.yHat{subjInd}(:,tInd),1);
        resSubj.auc(subjInd,:,tInd) = auc(1);
        resSubj.auc_CI5(subjInd,:,tInd) = nan;
        resSubj.auc_CI95(subjInd,:,tInd) = nan;
        %     [~,~,~,auc] = perfcurve(resSubj.y{subjInd},resSubj.yHat{subjInd},1,'NBOOT',2^10);
        %     resSubj.auc(subjInd) = auc(1);
        %     resSubj.auc_CI5(subjInd) = auc(2);
        %     resSubj.auc_CI95(subjInd) = auc(3);
        
        [~,P,~,STATS] = ttest(resSubj.yHat{subjInd}(resSubj.y{subjInd}==1,tInd),resSubj.yHat{subjInd}(resSubj.y{subjInd}==2,tInd));
        resSubj.distT(subjInd,:,tInd) = STATS.tstat;
        resSubj.distT_p(subjInd,:,tInd) = P;
        
%         resSubj.distMean(subjInd,:,tInd) = mean(resSess.distMean(subjInd,:,tInd),2);
%         resSubj.distStd(subjInd,:,tInd) = std(resSess.distMean(subjInd,:,tInd),[],2);
    end
    resSubj.distMean(subjInd,:,:) = mean(resSess.distMean(subjInd,:,:),2);
    resSubj.distStd(subjInd,:,:) = std(resSess.distMean(subjInd,:,:),[],2);
end
resSubj.nVoxOrig = round(median(resSess.nVoxOrig,2));
resSubj.nVox = round(median(resSess.nVox,2));
resSubj.svmSpace = resSess.svmSpace(:,1);

resGroup.y = cell2mat(resSubj.y);
resGroup.yHat = cell2mat(resSubj.yHat);
resGroup.nObs = sum(resSubj.nObs,1);
if ~isempty(resSubj.hit)
    resGroup.hit = sum(resSubj.hit,1);
    
    resGroup.acc = resGroup.hit./resGroup.nObs;
    [~,pci] = binofit(resGroup.hit,resGroup.nObs,0.1);
    resGroup.acc_CI5 = pci(:,1);
    resGroup.acc_CI95 = pci(:,2);
    [~,pci] = binofit(resGroup.nObs/2,resGroup.nObs,0.1);
    resGroup.acc_thresh = pci(:,2);
    resGroup.acc_p = binocdf(resGroup.hit,resGroup.nObs,0.5,'upper');
    [~,P,~,STATS] = ttest(resSubj.acc,0.5,'tail','right');
    resGroup.acc_T = STATS.tstat;
    resGroup.acc_P = P;
    [P,~,STATS] = signrank(resSubj.acc,0.5,'tail','right');
    resGroup.acc_wilcoxonSignedrank = STATS.signedrank;
    resGroup.acc_wilcoxonP = P;
else
    resGroup.hit = [];
    
    resGroup.acc = [];
    resGroup.acc_CI5 = [];
    resGroup.acc_CI95 = [];
    resGroup.acc_thresh = [];
    resGroup.acc_p = [];
    resGroup.acc_T = [];
    resGroup.acc_P = [];
    resGroup.acc_wilcoxonSignedrank = [];
    resGroup.acc_wilcoxonP = [];
end

resGroup.auc = nan(1,1,size(resGroup.yHat,2));
resGroup.auc_CI5 = nan(1,1,size(resGroup.yHat,2));
resGroup.auc_CI95 = nan(1,1,size(resGroup.yHat,2));
resGroup.auc_T = nan(1,1,size(resGroup.yHat,2));
resGroup.auc_P = nan(1,1,size(resGroup.yHat,2));
resGroup.auc_wilcoxonSignedrank = nan(1,1,size(resGroup.yHat,2));
resGroup.auc_wilcoxonP = nan(1,1,size(resGroup.yHat,2));
resGroup.distT = nan(1,1,size(resGroup.yHat,2));
resGroup.distT_p = nan(1,1,size(resGroup.yHat,2));
resGroup.distT_T = nan(1,1,size(resGroup.yHat,2));
resGroup.distT_P = nan(1,1,size(resGroup.yHat,2));
resGroup.distMean_T = nan(1,1,size(resGroup.yHat,2));
resGroup.distMean_P = nan(1,1,size(resGroup.yHat,2));
for tInd = 1:size(resGroup.yHat,2)
    [~,~,~,auc] = perfcurve(resGroup.y,resGroup.yHat(:,tInd),1,'NBOOT',2^10);
    resGroup.auc(tInd) = auc(1);
    resGroup.auc_CI5(tInd) = auc(2);
    resGroup.auc_CI95(tInd) = auc(3);
    [~,P,~,STATS] = ttest(resSubj.auc(:,:,tInd),0.5,'tail','right');
    resGroup.auc_T(tInd) = STATS.tstat;
    resGroup.auc_P(tInd) = P;
    [P,~,STATS] = signrank(resSubj.auc(:,:,tInd),0.5,'tail','right');
    resGroup.auc_wilcoxonSignedrank(tInd) = STATS.signedrank;
    resGroup.auc_wilcoxonP(tInd) = P;

    [~,P,~,STATS] = ttest(resGroup.yHat(resGroup.y==1,tInd),resGroup.yHat(resGroup.y==2,tInd),'tail','right');
    resGroup.distT(tInd) = STATS.tstat;
    resGroup.distT_p(tInd) = P;
    [~,P,~,STATS] = ttest(resSubj.distT(:,:,tInd),0,'tail','right');
    resGroup.distT_T(tInd) = STATS.tstat;
    resGroup.distT_P(tInd) = P;
    [~,P,~,STATS] = ttest(resSubj.distMean(:,:,tInd),0,'tail','right');
    resGroup.distMean_T(tInd) = STATS.tstat;
    resGroup.distMean_P(tInd) = P;
end
resGroup.nVoxOrig = round(median(resSubj.nVoxOrig,1));
resGroup.nVox = round(median(resSubj.nVox,1));
resGroup.svmSpace = resSubj.svmSpace{1};
