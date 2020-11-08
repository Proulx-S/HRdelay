function [hitRate, svmStruct,sortInd,w,A,mi,data,xValCorr] = runSVM(d,p,verbose,sortIndFix)

data = [];
global err
svmStruct = cell(1,1);
if ~exist('verbose','var')
    verbose = 0;
end
if verbose
    display(['Running svm at ' num2str(length(p.featLevList)) ' feature selection levels'])
end
if ~isfield(p,'rotTo')
    p.rotTo = pi/4;
end

if ischar(p.k)
    switch p.k
        case {'auto','autoCmplt'}
            k = length(d.crossVal)/2;
        case 'autoRun'
            k = length(d.crossVal)/2/str2double(p.split(1));
        otherwise
        error('');
    end
else
    k = p.k;
    p.k = 'randK';
end

%Loop over cross-validation folds
dOrig = d;
hitRate = nan(k,length(p.featLevList),length(p.featLevList));

w = nan(length(p.featLevList),length(p.featLevList),p.featLevList(end),k);
A = nan(length(p.featLevList),length(p.featLevList),p.featLevList(end),k);
% data = nan(length(p.featLevList),p.featLevList(end),k);
sortInd = nan(k,p.nFeatures);
mi = nan(k,p.nFeatures);
xValCorr.Rsquared = nan(k,length(p.featLevList));
xValCorr.RHO = nan(k,length(p.featLevList));

%% Loop over cross-validation folds
for fold = 1:k
    if verbose
        tic
    end
    if verbose
        display(['Fold ' num2str(fold) '/' num2str(k)])
    end
    % Initiate some stuff
    d = dOrig;
    
    
%     switch p.preSelection
%         case 'allSess'
%             error('double-check that')
%         case 'eqSess'
% %             error('double-check that')
% %             nRuns = size(d.xData,1);
% %             % Preselect data (randomly remove runs such that each session has equal number of runs)
% %             tmp = unique(d.sessionLabel);
% %             for i = 1:length(tmp)
% %                 sessionN(i) = length(find(d.sessionLabel==i));
% %             end
% %             run2remove = (sessionN-min(sessionN))./2;
% %             sessionLabel = d.sessionLabel(1:end/2);
% %             ind2rm = [];
% %             for i = 1:length(tmp)
% %                 ind2rm_tmp = find(sessionLabel==i);
% %                 ind2rm_tmp = ind2rm_tmp(randperm(length(ind2rm_tmp),run2remove(i)));
% %                 ind2rm = [ind2rm ind2rm_tmp];
% %             end
% %             d.xData([ind2rm; nRuns/2+ind2rm],:) = [];
% %             d.normData(ind2rm,:) = [];
% %             d.label([ind2rm; nRuns/2+ind2rm]) = [];
% %             d.sessionLabel([ind2rm; nRuns/2+ind2rm]) = [];
% %             d.crossVal([ind2rm; nRuns/2+ind2rm]) = [];
%         case 'sess1'
%             nRuns = size(d.xData,1);
%             ind2rm = find(d.sessionLabel(1:end/2)~=1);
%             d.xData([ind2rm; nRuns/2+ind2rm],:) = [];
%             d.normData(ind2rm,:) = [];
%             d.label([ind2rm; nRuns/2+ind2rm]) = [];
%             d.sessionLabel([ind2rm; nRuns/2+ind2rm]) = [];
%             d.crossVal([ind2rm; nRuns/2+ind2rm]) = [];
%         case 'sess2'
%             nRuns = size(d.xData,1);
%             ind2rm = find(d.sessionLabel(1:end/2)~=2);
%             d.xData([ind2rm; nRuns/2+ind2rm],:) = [];
%             d.normData(ind2rm,:) = [];
%             d.label([ind2rm; nRuns/2+ind2rm]) = [];
%             d.sessionLabel([ind2rm; nRuns/2+ind2rm]) = [];
%             d.crossVal([ind2rm; nRuns/2+ind2rm]) = [];
%     end
    
    
    
    %     curFoldLabel = foldsList(fold);
    wFold = w(:,:,:,fold);
    Afold = A(:,:,:,fold);
    %     dataFold = data(:,:,fold);
    
    
    % Break complex data down to amp and delay
    d2.delay(:,:,1) = angle(d.xData(1:end/2,:));
    d2.delay(:,:,2) = angle(d.xData(end/2+1:end,:));
    d2.delay(:,:,3) = angle(d.normData);
    d2.amp(:,:,1) = abs(d.xData(1:end/2,:));
    d2.amp(:,:,2) = abs(d.xData(end/2+1:end,:));
    d2.amp(:,:,3) = abs(d.normData);
    
    switch p.fixDelay
        case 'vox'
            d2.delay(:,:,1) = wrapToPi(d2.delay(:,:,1)-repmat(circ_mean(d2.delay(:,:,3),[],1),[size(d2.delay,1) 1]));
            d2.delay(:,:,2) = wrapToPi(d2.delay(:,:,2)-repmat(circ_mean(d2.delay(:,:,3),[],1),[size(d2.delay,1) 1]));
            d2.delay(:,:,3) = wrapToPi(d2.delay(:,:,3)-repmat(circ_mean(d2.delay(:,:,3),[],1),[size(d2.delay,1) 1]));
            [d2.amp,~] = pol2cart(d2.delay,d2.amp);
        case 'pop'
            error('not implemented yet')
        case 'free'
        otherwise
            error('specificy that')
    end
    d2.sessionLabel(:,:,1) = d.sessionLabel(1:end/2,:);
    d2.sessionLabel(:,:,2) = d.sessionLabel(end/2+1:end,:);
    d2.sessionLabel(:,:,3) = d.sessionLabel(1:end/2,:);
    d2.label(:,:,1) = d.label(1:end/2,:);
    d2.label(:,:,2) = d.label(end/2+1:end,:);
    d2.label(:,:,3) = ones(size(d.label(1:end/2,:))).*3;
    if ~strcmp(p.k,'randK')
        d2.crossVal(:,:,1) = d.crossVal(1:end/2,:);
        d2.crossVal(:,:,2) = d.crossVal(end/2+1:end,:);
        d2.crossVal(:,:,3) = nan(size(d.crossVal(1:end/2,:)));
    end
    
    
    %     close all
    %     figure('windowStyle','docked');
    %     tmp = d2.amp;
    %     tmp = reshape(tmp,[size(tmp,1)*size(tmp,3) size(tmp,2)])';
    %     tmp = corr(tmp);
    %     imagesc(tmp,[0 1])
    %
    %
    %
    %     close all
    %     figure('windowStyle','docked');
    %     tmp = abs(cat(1,d.xData,d.normData))';
    %     tmp = corr(tmp);
    %     imagesc(tmp,[0 1])
    
    
    %     d.amp = abs(d.xData);
    %     d.delay = angle(d.xData);
    
    
    
    
    %     % Normalization to plaid
    %     if p.normToPlaid
    %         %Amp
    %         d.amp = d.amp - repmat(abs(d.normData),[2 1]);
    %         %Delay
    %         d.delay = wrapToPi(d.delay - repmat(angle(d.normData),[2 1]));
    %     end
    
    
    
    %     % Define train and test data (all sample gets to be in the testing set once; not finished)
    %     labels = p2.sessionLabel;
    %     sessionLabels = p2.sessionLabel(:,:,1);
    %     sessions = unique(sessionLabels);
    %
    %     teInd = {};
    %     while ~all(isnan(sessionLabels))
    %         sessions = unique(sessionLabels);
    %         sessions = sessions(~isnan(sessions));
    %         for i = 1:length(sessions)
    %             curInd = find(sessionLabels(:,:,1)==sessions(i));
    %             curTeInd(i) = curInd(randperm(length(curInd),1));
    %         end
    %         teInd{end+1} = curTeInd;
    %         sessionLabels(curTeInd) = nan;
    %         curTeInd = [];
    %     end
    %     teInd{1}
    %     teInd{2}
    %     teInd{3}
    
    switch p.k
        case {'auto','autoCmplt'}
            error('teInd and trInd badly defined')
            teInd = squeeze(d2.crossVal(fold,1,1:2));
            d2tr = d2;
            allFields = fields(d2tr);
            for i = 1:length(allFields)
                d2te.(allFields{i}) = d2tr.(allFields{i})(teInd,:,:);
                d2tr.(allFields{i})(teInd,:,:) = [];
            end
        case 'autoRun'
            teInd1 = d2.crossVal(:,1,1)==fold;
            teInd2 = d2.crossVal(:,1,2)==fold;
            allFields = fields(d2);
            for i = 1:length(allFields)
                d2te.(allFields{i})(:,:,1) = d2.(allFields{i})(teInd1,:,1);
                d2te.(allFields{i})(:,:,2) = d2.(allFields{i})(teInd2,:,2);
                tmp1 = d2.(allFields{i})(~teInd1,:,1);
                tmp2 = d2.(allFields{i})(~teInd2,:,2);
                d2tr.(allFields{i}) = cat(3,tmp1,tmp2); clear tmp1 tmp2
            end
        case 'sessBal'
            error('double-check that')
        case 'randK'
            % Define train and test data (purely random selection of one sample per session for the testing set, puts more weights on sessions with fewer sample)
            sessions = unique(d2.sessionLabel(:,:,1));
            for i = 1:length(sessions)
                curInd = find(d2.sessionLabel(:,:,1)==sessions(i));
                teInd(i) = curInd(randperm(length(curInd),1));
            end
            d2tr = d2;
            allFields = fields(d2tr);
            for i = 1:length(allFields)
                d2te.(allFields{i}) = d2tr.(allFields{i})(teInd,:,:);
                d2tr.(allFields{i})(teInd,:,:) = [];
            end
        otherwise
            error('double-check that')
    end
    
    
    
    
    
    %     [tr.amp,te.amp,~,~,~,~] = splitTrainAndTest(d,curFoldLabel,'amp');
    %     [tr.delay,te.delay,tr.label,tr.sessionLabel,te.label,te.sessionLabel,tr.label1,tr.label2] = splitTrainAndTest(d,curFoldLabel,'delay');
    
    
    
    %     % Regress out session effect
    %     close all
    %     figure('windowStyle','docked');
    %     tmp = d2.amp;
    %     tmp = cat(1,tmp(:,:,1),tmp(:,:,2),tmp(:,:,3));
    %     tmp = corr(tmp');
    %     imagesc(tmp,[0 1]); colorbar
    %     title('before')
    %     saveas(gca,fullfile(p.dataDirOut,[p.subj '_corrMat_bf.fig']))
    %     saveas(gca,fullfile(p.dataDirOut,[p.subj '_corrMat_bf.jpg']))
    %
    %     [tmp,X,curX] = regressOutSessionVoxWise2(d2.amp,d2.sessionLabel,1);
    %
    %     figure('windowStyle','docked');
    %     tmp = cat(1,tmp(:,:,1),tmp(:,:,2),tmp(:,:,3));
    %     tmp = corr(tmp');
    %     imagesc(tmp,[0 1]); colorbar
    %     title('after')
    %     saveas(gca,fullfile(p.dataDirOut,[p.subj '_corrMat_af.fig']))
    %     saveas(gca,fullfile(p.dataDirOut,[p.subj '_corrMat_af.jpg']))
    %
    %     figure('windowStyle','docked');
    %     imagesc(X); colormap gray
    %     saveas(gca,fullfile(p.dataDirOut,[p.subj '_designMat.fig']))
    %     saveas(gca,fullfile(p.dataDirOut,[p.subj '_designMat.jpg']))
    %     figure('windowStyle','docked');
    %     imagesc(curX); colormap gray
    %     saveas(gca,fullfile(p.dataDirOut,[p.subj '_nuisanceMat.fig']))
    %     saveas(gca,fullfile(p.dataDirOut,[p.subj '_nuisanceMat.jpg']))
    
    
    
    if p.regSession
        error('double-check that')
        %estimate session effect from training set
        [betaHat,X] = estimateSessionEffect(d2tr.amp,d2tr.sessionLabel);
        firsSessInd = size(X,2)-(length(unique(d2tr.sessionLabel))-1)+1;
        sessionEffect = betaHat(:,firsSessInd:end)';
        
        %remove session effect
        for sess = 2:length(unique(d2tr.sessionLabel))
            %from training set
            sessInd = d2tr.sessionLabel(:,:,1)==sess;
            d2tr.amp(sessInd,:,:) = d2tr.amp(sessInd,:,:)-repmat(sessionEffect,[length(find(sessInd)) 1 3]);
            %from testing set
            sessInd = d2te.sessionLabel(:,:,1)==sess;
            d2te.amp(sessInd,:,:) = d2te.amp(sessInd,:,:)-repmat(sessionEffect,[length(find(sessInd)) 1 3]);
        end
    end
    
    
    
    % Normalize (previously done with newNormalization.m)
    allFields = fields(d2tr);
    for i = 1:length(allFields)
        tr.(allFields{i}) = cat(1,d2tr.(allFields{i})(:,:,1),d2tr.(allFields{i})(:,:,2));
    end
    allFields = fields(d2te);
    for i = 1:length(allFields)
        te.(allFields{i}) = cat(1,d2te.(allFields{i})(:,:,1),d2te.(allFields{i})(:,:,2));
    end
    
    trPreNorm = tr;
    tePreNorm = te;
    %Delay
    %remove mean phase and rotate
    phase_shift = circ_mean(tr.delay,[],1)-p.rotTo;
    tr.delay = wrapToPi(tr.delay-repmat(phase_shift,size(tr.delay,1),1));
    te.delay = wrapToPi(te.delay-repmat(phase_shift,size(te.delay,1),1));
    %Amp
    %shift
    amp_shift = mean(tr.amp,1);
    tr.amp = tr.amp-repmat(amp_shift,size(tr.amp,1),1);
    te.amp = te.amp-repmat(amp_shift,size(te.amp,1),1);
    %scale
    amp_scale = std(tr.amp,[],1);
    tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
    te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
    
    %Sorting for feature slection (previously done with newFeatSorting.m)
    %compute sorting
    if ~exist('sortIndFix','var')
        switch p.keepInfo
            case 'm'
                p.sort.feat = computeMI(tr.amp,tr.label==1,tr.label==2)';
                %                 [~,p.sort.indFeat] = sort(p.sort.feat,'descend');
            case 'p'
                error('not implemented yet')
        end
    else
        p.sort.feat = nan(1,p.nFeatures);
        %         p.sort.indFeat = sortIndFix;
    end
    
    %Sorting for functional ROI
    p.sort.fun = p.funcROI.vec.stats;
    [~,p.sort.indFun] = sort(p.sort.fun,'descend');
    
    
    %     %% Look at session effect
    %     keyboard
    %     tmp = tr.amp(:,p.funcROI.vec.sortInd);
    %     tmp1 = tmp(1:end/2,:);
    %     tmp2 = tmp(end/2+1:end,:);
    %     ind = mean(tmp2-tmp1,1)>0;
    %     tmp(1:end/2,ind) = tmp1(:,ind);
    %     tmp(end/2+1:end,ind) = tmp2(:,ind);
    %     tmp(1:end/2,~ind) = tmp2(:,~ind);
    %     tmp(end/2+1:end,~ind) = tmp1(:,~ind);
    %
    %     imagesc(tmp)
    
    
    trPreFunSort = tr;
    tePreFunSort = te;
    
    
    
    
    %% Loop over functional ROI levels
    for funLevInd = 1:length(p.featLevList)
        %% Loop over feature selection levels
        %apply fun sorting (and select fun)
        trFeatSort.amp = trPreFunSort.amp(:,p.sort.indFun((1:p.featLevList(funLevInd))));
        trFeatSort.delay = trPreFunSort.delay(:,p.sort.indFun((1:p.featLevList(funLevInd))));
        teFeatSort.amp = tePreFunSort.amp(:,p.sort.indFun((1:p.featLevList(funLevInd))));
        teFeatSort.delay = tePreFunSort.delay(:,p.sort.indFun((1:p.featLevList(funLevInd))));
        curFun = p.sort.fun(p.sort.indFun((1:p.featLevList(funLevInd))));
        curFeat = p.sort.feat(p.sort.indFun((1:p.featLevList(funLevInd))));
        
        %sort feat
        [~,curIndFeat] = sort(curFeat,'descend');
        
        %apply feat sorting
        trFeatSort.amp = trFeatSort.amp(:,curIndFeat);
        trFeatSort.delay = trFeatSort.delay(:,curIndFeat);
        teFeatSort.amp = teFeatSort.amp(:,curIndFeat);
        teFeatSort.delay = teFeatSort.delay(:,curIndFeat);
        curFun = curFun(curIndFeat);
        curFeat = curFeat(curIndFeat);
        
        for featLevInd = 1:length(p.featLevList)
            if p.featLevList(featLevInd)>length(curFun)
                break
            end
            %select feat
            tr.amp = trFeatSort.amp(:,1:p.featLevList(featLevInd));
            tr.delay = trFeatSort.delay(:,1:p.featLevList(featLevInd));
            te.amp = teFeatSort.amp(:,1:p.featLevList(featLevInd));
            te.delay = teFeatSort.delay(:,1:p.featLevList(featLevInd));
            curFun(1:p.featLevList(featLevInd));
            curFeat(1:p.featLevList(featLevInd));
            
            
            % Cross-validated pattern correlation
            if p.doCrossVal
                switch p.keepInfo
                    case 'm'
                        trPat = mean(tr.amp(1:end/2,:) - tr.amp(end/2+1:end,:),1);
                        tePat = mean(te.amp(1:end/2,:) - te.amp(end/2+1:end,:),1);
                        xValCorr.RHO(fold,featLevInd) = corr(trPat',tePat');
                        xValCorr.Rsquared(fold,featLevInd) = 100*(1 - sum((tePat-trPat).^2) / sum((tePat-mean(tePat,2)).^2));
                    case 'p'
                        trPat = circ_mean(wrapToPi(tr.delay(1:end/2,:) - tr.delay(end/2+1:end,:)),[],1);
                        tePat = circ_mean(wrapToPi(te.delay(1:end/2,:) - te.delay(end/2+1:end,:)),[],1);
                        xValCorr.RHO(fold,featLevInd) = circ_corrcc(trPat', tePat');
                        xValCorr.Rsquared(fold,featLevInd) = 100*(1 - sum(wrapToPi(tePat-trPat).^2) / sum(wrapToPi(tePat-circ_mean(tePat,[],2)).^2));
                end
            end
            
            
            % SVM
            if p.doSVM
                %Train and test
                switch p.keepInfo
                    case 'm'
                        complexFlag = 0;
                        trFinal.data = tr.amp;
                        trFinal.label = tr.label;
                        teFinal.data = te.amp;
                        teFinal.label = te.label;
                    case 'p'
                        complexFlag = 1;
                        [X,Y] = pol2cart(tr.delay,1);
                        trFinal.data = [X Y];
                        trFinal.label = tr.label;
                        [X,Y] = pol2cart(te.delay,1);
                        teFinal.data = [X Y];
                        teFinal.label = te.label;
                end
                
                
%                 %# construct predicting attributes and target class
%                 x = trFinal.data;
%                 y = cellstr(num2str(trFinal.label));
% %                 y = (trFinal.label-1)*2-1;
%                 
%                 %# train classification decision tree
% %                 t = classregtree(x, y, 'method','classification','prune','off');
%                 cost = diag([nan nan]);
%                 cost(~isnan(cost)) = 0.5;
%                 cost(isnan(cost)) = 0;
%                 
%                 t = classregtree(x, y,'minparent',1,'cost',cost)
%                 t = classregtree(x, y,'minparent',1,'cost',cost,'splitcriterion','deviance')
%                 view(t)
%                 
%                 %# test
%                 yPredicted = eval(t, x);
%                 cm = confusionmat(y,yPredicted);           %# confusion matrix
%                 N = sum(cm(:));
%                 err = ( N-sum(diag(cm)) ) / N;             %# testing error
%                 
% %                 %# prune tree to avoid overfitting
% %                 tt = prune(t, 'level',3);
% %                 view(tt)
%                 
%                 %# predict a new unseen instance
%                 xte = teFinal.data;
%                 yte = cellstr(num2str(teFinal.label));
%                 prediction = eval(t, xte);    %# pred = 'Japan'
                
                
                [hitRate(fold,funLevInd,featLevInd),svmStruct] = trainNtestSVM(trFinal.data,trFinal.label,teFinal.data,teFinal.label,p);
                %             [hitRate(fold,featLevInd),svmStruct] = trainNtestSVM(trFinal.data,trFinal.label,teFinal.data,teFinal.label,p);
                
                %Extract w
                alpha = zeros(length(svmStruct.GroupNames),1);
                alpha(svmStruct.SupportVectorIndices,1) = svmStruct.Alpha;
                y = (svmStruct.GroupNames-1)*2-1;
                x = trFinal.data;
                wTmp = zeros(1,size(x,2));
                for i = 1:size(x,1)
                    wTmp = wTmp + alpha(i)*y(i)*x(i,:);
                end
                
                %Extract s
                s = wTmp*x';
                
                %Extract A
                Atmp = wTmp*cov(x)/cov(s);
                
                %Compile
                if complexFlag
                    wTmp = complex(wTmp(1:end/2),wTmp(end/2+1:end));
                    Atmp = complex(Atmp(1:end/2),Atmp(end/2+1:end));
                end
                wFold(funLevInd,featLevInd,1:p.featLevList(featLevInd)) = wTmp;
                Afold(funLevInd,featLevInd,1:p.featLevList(featLevInd)) = Atmp;
                %             dataFold(featLevInd,1:p.featLevList(featLevInd)) = mean(tr_origSorted_cur,1);
            end
        end
    end
%     if p.doSVM
%         %Sort back weights to original and store
%         w(:,:,fold) = sort_back(wFold,repmat(sortInd(fold,:),size(wFold,1),1),2);
%         A(:,:,fold) = sort_back(Afold,repmat(sortInd(fold,:),size(Afold,1),1),2);
%         %         data(:,:,fold) = sort_back(dataFold,repmat(sortInd(fold,:),size(dataFold,1),1),2);
%     end
    if verbose
        toc
    end
end

%Average over cross-validation folds
hitRate = mean(hitRate,1);
xValCorr.Rsquared = mean(xValCorr.Rsquared,1);
xValCorr.RHO = mean(xValCorr.RHO,1);


