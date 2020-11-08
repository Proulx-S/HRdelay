function [hitRate,allInd2keep] = runSVM_RFE(d,p)

%Train and test SVM at each cross-validation fold
tmpFeatLevList = sort(p.featLevList,'descend');
foldsList = unique(d.crossVal);
hitRate = nan(length(foldsList),length(tmpFeatLevList));
for fold = 1:length(foldsList)
    display(['Fold ' num2str(fold)])
    curFoldLabel = foldsList(fold);
    
    %Define train and test data
    trainInd = true(size(d.xData,1),1);
    trainInd(d.crossVal==curFoldLabel) = false;
    testInd = false(size(d.xData,1),1);
    testInd(d.crossVal==curFoldLabel) = true;
    tr = d.xData(trainInd,:);
    trLabel = d.label(trainInd);
    te = d.xData(testInd,:);
    teLabel = d.label(testInd);

    %Tranform data appropriately
    selectedFeatures = 1:p.nFeatures;
    [te,tr] = transformSelectedFeatures(te,tr,selectedFeatures,p);
        
    %Run RFE
    for curFeatLev = 1:length(tmpFeatLevList);
        if tmpFeatLevList(curFeatLev)~=p.nFeatures
            %Get low weigthed features
            [ind2keep] = eliminateFeature(tr,trLabel,p,tmpFeatLevList(curFeatLev));
            %Store features
            allInd2keep{fold,curFeatLev} = ind2keep;
            %Rearange for complex values
            switch p.learningVariable % 'r' 'i' 'ri' 'm' 'p' 'mp'
                case 'ri'
                    ind2keep = [ind2keep ind2keep+tmpFeatLevList(curFeatLev)];
                case {'r','i'}
                case {'m','p'}
                case {'mp'}
                    error('not implemented')
                otherwise
                    error('p.curParam badly specifie')
            end
            %Eliminate low weigthed features
            tr = tr(:,ind2keep);
            te = te(:,ind2keep);
        end
        %Train and test on reduced data set
        hitRate(fold,curFeatLev) = trainNtestSVM(tr,trLabel,te,teLabel,p);
    end
end
%Average across folds
hitRate = mean(hitRate,1);
%Flip results since we did it the reverse as asked in p.featLevlist
if any(tmpFeatLevList~=p.featLevList)
    hitRate = fliplr(hitRate);
end
