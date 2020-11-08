function [hitRate, svmStruct,sortInd] = runSVMcomplex(d,p,verbose,sortIndFix)
sortInd = [];
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
    [trAll,teAll] = newNormalization(tr,te,p);
    %Sorting for feature slection
    [trAll,teAll,sortInd(fold,:)] = newFeatSorting(trAll,teAll,trLabel1,trLabel2,p);
        
    for featLevInd = 1:length(p.featLevList)
        %Feature selection
        tr = trAll(:,1:p.featLevList(featLevInd));
        te = teAll(:,1:p.featLevList(featLevInd));
        %Transform for complex svm
        tr = [real(tr) imag(tr)];
        te = [real(te) imag(te)];
%         tr = [real(tr)];
%         te = [real(te)];
        %Train and test
        [hitRate(fold,featLevInd),svmStruct] = trainNtestSVM(tr,trLabel,te,teLabel,p);
%         toc
    end
end
%Average over cross-validation folds
hitRate = mean(hitRate,1);



