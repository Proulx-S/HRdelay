function [hitRate, svmStruct,sortInd,w,A,mi,data,Rsquared] = runSVM(d,p,verbose,sortIndFix)

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

w = nan(length(p.featLevList),p.featLevList(end),length(foldsList));
A = nan(length(p.featLevList),p.featLevList(end),length(foldsList));
data = nan(length(p.featLevList),p.featLevList(end),length(foldsList));
sortInd = nan(length(foldsList),p.nFeatures);
mi = nan(length(foldsList),p.nFeatures);
Rsquared = nan(length(foldsList),length(p.featLevList));
for fold = 1:length(foldsList)
    d = dOrig;
    curFoldLabel = foldsList(fold);
    wFold = w(:,:,fold);
    Afold = A(:,:,fold);
    dataFold = data(:,:,fold);
    
    
    
    %Normalization to plaid
    if p.normToPlaid
        if p.filterData
            switch p.keepInfo
                case 'p'
                    phase = wrapToPi(angle(d.xData) - repmat(angle(d.normData),[2 1]));
                    [X,Y] = pol2cart(phase,1);
                    d.xData = complex(X,Y); clear X Y
                case 'm'
                    d.xData = abs(d.xData) - repmat(abs(d.normData),[2 1]);
                otherwise
                    error('not implemented for that')
            end
        else
            error('not sure who to do that yet')
        end
    end
    
    %Define train and test data
    [tr_orig,te_orig,trLabel,teLabel,trLabel1,trLabel2] = splitTrainAndTest(d,curFoldLabel);
    
    %% Loop over feature selection levels
    if verbose
        display(['Fold ' num2str(fold) '/' num2str(p.k)])
    end
    
    %Normalize
    [trAll,teAll] = newNormalization(tr_orig,te_orig,p); % Here we filter according to p.keepInfo and normalize
    
    %Sorting for feature slection
    if ~exist('sortIndFix','var')
        [trAll,teAll,sortInd(fold,:),mi(fold,:)] = newFeatSorting(trAll,teAll,trLabel1,trLabel2,p);
    else
        if length(find(size((sortIndFix))>1))>1
            sortInd(fold,:) = sortIndFix(fold,:);
        else
            sortInd(fold,:) = sortIndFix;
        end
        trAll = trAll(:,sortInd(fold,:));
        teAll = teAll(:,sortInd(fold,:));
        mi = [];
    end
    tr_origSorted = tr_orig(:,sortInd(fold,:));  % Just to be able to output delay and amp in a way that matches svm
    
    for featLevInd = 1:length(p.featLevList)
%         if featLevInd == length(p.featLevList)
%             keyboard
%         end
        %Feature selection
        tr = trAll(:,1:p.featLevList(featLevInd));
        te = teAll(:,1:p.featLevList(featLevInd));
        
        tr_origSorted_cur = tr_origSorted(:,1:p.featLevList(featLevInd));
        
        
        %Transform for complex svm
        if isreal(tr)
            complexFlag = 0;
        else
            complexFlag = 1;
            tr = [real(tr) imag(tr)];
            te = [real(te) imag(te)];
        end
%         if featLevInd==length(p.featLevList)
%             keyboard
%         end

        %Cross-validated pattern correlation
        if p.doCrossVal
            switch p.keepInfo
                case 'p'
                    [trP1,~] = cart2pol(tr(1:end/2,1:end/2),tr(1:end/2,end/2+1:end));
                    [trP2,~] = cart2pol(tr(end/2+1:end,1:end/2),tr(end/2+1:end,end/2+1:end));
                    trP = circ_mean(wrapToPi(trP2-trP1),[],1);
                    [teP1,~] = cart2pol(te(1:end/2,1:end/2),te(1:end/2,end/2+1:end));
                    [teP2,~] = cart2pol(te(end/2+1:end,1:end/2),te(end/2+1:end,end/2+1:end));
                    teP = wrapToPi(teP2-teP1);
                    Rsquared(fold,featLevInd) = 100*(1 - sum((wrapToPi(teP-trP)).^2) / sum((wrapToPi(teP-circ_mean(teP,[],2))).^2));
                case 'm'
                    trP = mean(tr(1:end/2,:) - tr(end/2+1:end,:),1);
                    teP = te(1:end/2,:) - te(end/2+1:end,:);
                    Rsquared(fold,featLevInd) = corr(trP',teP');
%                     Rsquared(fold,featLevInd) = 100*(1 - sum((teP-trP).^2) / sum((teP-mean(teP,2)).^2));
            end
        end
        
        
        %Train and test
        if p.doSVM
            [hitRate(fold,featLevInd),svmStruct] = trainNtestSVM(tr,trLabel,te,teLabel,p);
            %Extract w
            alpha = zeros(length(svmStruct.GroupNames),1);
            alpha(svmStruct.SupportVectorIndices,1) = svmStruct.Alpha;
            y = (svmStruct.GroupNames-1)*2-1;
            x = tr;
            wTmp = zeros(1,size(x,2));
            for i = 1:size(x,1)
                wTmp = wTmp + alpha(i)*y(i)*x(i,:);
            end
            
            %Extract s
            s = wTmp*x';
            
            %Extract A
            Atmp = wTmp*cov(x)/cov(s);
            
            if complexFlag
                wTmp = complex(wTmp(1:end/2),wTmp(end/2+1:end));
                Atmp = complex(Atmp(1:end/2),Atmp(end/2+1:end));
            end
            
            wFold(featLevInd,1:p.featLevList(featLevInd)) = wTmp;
            Afold(featLevInd,1:p.featLevList(featLevInd)) = Atmp;
            %         toc
            
            dataFold(featLevInd,1:p.featLevList(featLevInd)) = mean(tr_origSorted_cur,1);
        end
        
    end
    if p.doSVM 
        %Sort back weights to original and store
        w(:,:,fold) = sort_back(wFold,repmat(sortInd(fold,:),size(wFold,1),1),2);
        A(:,:,fold) = sort_back(Afold,repmat(sortInd(fold,:),size(Afold,1),1),2);
        data(:,:,fold) = sort_back(dataFold,repmat(sortInd(fold,:),size(dataFold,1),1),2);
    end
    
    %         toc
end

%Average over cross-validation folds
hitRate = mean(hitRate,1);
Rsquared = mean(Rsquared,1);


