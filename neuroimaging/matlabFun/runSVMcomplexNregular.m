function [hitRate, svmStruct,proportionPhase] = runSVM4(d,p,verbose,sortIndFix)
global err
svmStruct = cell(1,1);
if ~exist('verbose','var')
    verbose = 0;
end

if verbose
    display(['Running svm at ' num2str(length(p.featLevList)) ' feature selection levels'])
end

%Loop over cross-validation folds
foldsList = unique(d.crossVal);
dOrig = d;
hitRate = nan(length(foldsList),length(p.featLevList));
nPind = nan(length(foldsList),length(p.featLevList));
nMind = nan(length(foldsList),length(p.featLevList));

for fold = 1:length(foldsList)
    d = dOrig;
    curFoldLabel = foldsList(fold);
    
    %Define train and test data
    [tr,te,trLabel,teLabel,trLabel1,trLabel2] = splitTrainAndTest(d,curFoldLabel);

    %% Loop over feature selection levels
    if verbose
        display(['Fold ' num2str(fold) '/' num2str(p.k)])
    end
    
    %Normalize
    p.rotTo = 45;
    p.filterData = 1;
    p.keepInfo = 'p';
    [trP,teP] = newNormalization(tr,te,p);
    
    p.rotTo = 0;
    p.filterData = 1;
    p.keepInfo = 'm';
    [trM,teM] = newNormalization(tr,te,p);
    trM = real(trM); teM = real(teM);
    
    %Compute MI
    miP = computeMIcomplex(trP,trLabel1,trLabel2);
    miM = computeMIreal(trM,trLabel1,trLabel2);
    miAll = [miP miM];
    [~,sortInd] = sort(miAll,'descend');
        
    for featLevInd = 1:length(p.featLevList)
        %Feature selection
        curSortInd = sortInd(1:p.featLevList(featLevInd));
        nPind(fold,featLevInd) = length(find(curSortInd<=length(miP)));
        nMind(fold,featLevInd) = length(curSortInd>length(miP));
        trPcur = trP(:,curSortInd(curSortInd<=length(miP)));
        trMcur = trM(:,curSortInd(curSortInd>length(miP))-length(miP));
        tePcur = teP(:,curSortInd(curSortInd<=length(miP)));
        teMcur = teM(:,curSortInd(curSortInd>length(miP))-length(miP));
        
        
        %Transform for complex + regular svm
        trCur = cat(2,real(trPcur),imag(trPcur),trMcur);
        teCur = cat(2,real(tePcur),imag(tePcur),teMcur);
        %Train and test
        [hitRate(fold,featLevInd),svmStruct] = trainNtestSVM2(trCur,trLabel,teCur,teLabel,p);
%         toc
    end
end
%Average over cross-validation folds
hitRate = mean(hitRate,1);
proportionPhase = mean(nPind./(nMind+nPind),1);


