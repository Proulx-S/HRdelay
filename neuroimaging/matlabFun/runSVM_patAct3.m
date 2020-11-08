function out = runSVM_patAct3(d,p,verbose,rep,doPerm)

if ~exist('verbose','var')
    verbose = 0;
end
if ~exist('doPerm','var')
    doPerm = 0;
end

%% Permute data if needed
if doPerm
    d.xData = d.xData(randperm(size(d.xData,1)),:);
end

%% Define k-folding
if ischar(p.k)
    switch p.k
        %         case {'auto','autoCmplt'}
        %             k = length(d.crossVal)/2;
        case 'autoRun'
            d.crossVal = defineCrossVal(d,p);
            k = length(d.crossVal)/2/str2double(p.split(1));
        otherwise
            error('');
    end
else
    error('')
    k = p.k;
    p.k = 'randK';
end


%% Loop over cross-validation folds
dOrig = d;
% out.hitRate = nan(k,1);
% out.hitRate_d = nan(k,1);
% out.hitRate_pat = nan(k,1);
% out.hitRate_patRel = nan(k,1);
% 
% out.pat = nan(p.nObs,2);
% out.patRel = nan(p.nObs,2);

out.crossInfo = cell(length(p.trainInfo),length(p.testInfo));
out.crossInfo_dim1Train = p.trainInfo;
out.crossInfo_dim2Test = p.testInfo;
for trainInd = 1:length(p.trainInfo)
    for testInd = 1:length(p.testInfo)
%         display([num2str((trainInd-1)*length(p.testInfo)+testInd) '/' num2str(length(p.trainInfo)*length(p.testInfo))])
        hitRate = nan(k,1);
        decision = nan(p.nObs,1);
        for fold = 1:k
%             [out.crossInfo{trainInd,testInd}.hitRate(fold,1),out.crossInfo{trainInd,testInd}.d(:,fold)] = trainNtest5(dOrig,p,fold,k,rep,verbose,p.trainInfo{trainInd},p.testInfo{testInd});
            [hitRate(fold,1),curDecision,teInd] = trainNtest5(dOrig,p,fold,k,rep,verbose,p.trainInfo{trainInd},p.testInfo{testInd});
            decision(teInd) = curDecision;
        end
        out.crossInfo{trainInd,testInd}.hitRate = mean(hitRate,1);
        out.crossInfo{trainInd,testInd}.d = mean(decision,2);
    end
end

%Recompile
hitRate = nan(size(out.crossInfo));
dVal = nan([size(out.crossInfo) p.nObs]);
for trainInd = 1:length(p.trainInfo)
    for testInd = 1:length(p.testInfo)
        hitRate(trainInd,testInd) = out.crossInfo{trainInd,testInd}.hitRate;
        dVal(trainInd,testInd,:) = out.crossInfo{trainInd,testInd}.d;
    end
end
out.crossInfo = [];
out.crossInfo.hitRate = hitRate;
out.crossInfo.d = dVal;



