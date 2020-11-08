function [h,h2] = plotSVM(svm,h,testParam,testParamInd,h2)
if isfield(svm,'plotIt') && isnumeric(svm.plotIt) && svm.plotIt
    plotIt = svm.plotIt;
else
    plotIt = 1;
end

%% Plot SVM
if ~isfield(svm.p,'doSVM') || svm.p.doSVM
    %% Initiate plot
    if ~exist('h','var') || isempty(h)
        if plotIt
            h = figure('windowStyle','docked');
        else
            h = figure('windowStyle','docked','visible','off');
        end
    else
        if plotIt
            figure(h); hold off;
        else
            figure(h,'visible','off'); hold off;
        end
        
    end
    
    oldVersion = 0;
    if ~isfield(svm.p,'featLevList')
        svm.p.featLevList = 1:size(svm.r.hits,2);
    end
    if ~isfield(svm.p.binom,'hitRate')
        oldVersion = 1;
        svm.p.binom.hitRate = svm.p.binom.hits;
    end
    if ~isfield(svm.r,'hitRate')
        svm.r.hitRate = svm.r.hits;
    end
    
    includeNonParamDist = 0;
    if isfield(svm,'rP')
        includeNonParamDist = 1;
    end
    
    
    if oldVersion
        plot(svm.p.featLevList,svm.p.nObs/2*ones(size(svm.p.featLevList)),'-k'); hold on;
    else
        plot(svm.p.featLevList,0.5*ones(size(svm.p.featLevList)),'-k'); hold on;
    end
    
    %% Plot non-parametric distribution
    if includeNonParamDist
        if ~isfield(svm.pP,'featLevList')
            svm.pP.featLevList = 1:size(svm.rP.hits,2);
        end
        if ~isfield(svm.pP.binom,'hitRate')
            oldVersion = 1;
            svm.pP.binom.hitRate = svm.pP.binom.hits;
        end
        if ~isfield(svm.rP,'hitRate')
            svm.rP.hitRate = svm.rP.hits;
        end
        
        for i = 1:size(svm.rP.thresh,1)
            if ((size(svm.rP.thresh(i,:),2)==1) || (size(svm.rP.thresh(i,:),2)~=size(svm.r.hitRate,2)) || any(any(isnan(svm.rP.thresh))))
                plot(svm.pP.featLevList,svm.rP.thresh(i,:),'*r')
            else
                plot(svm.pP.featLevList,svm.rP.thresh(i,:),'--r')
            end
        end
        if ((size(svm.rP.thresh(i,:),2)==1) || (size(svm.rP.thresh(i,:),2)~=size(svm.r.hitRate,2)) || any(any(isnan(svm.rP.thresh))))
            plot(svm.pP.featLevList,nanmean(svm.rP.hitRate,1),'or')
        else
            plot(svm.pP.featLevList,mean(svm.rP.hitRate,1),'-r')
        end
        
        if ((size(svm.rP.thresh,2)==1) || (size(svm.rP.thresh,2)~=size(svm.r.hitRate,2)) || any(any(isnan(svm.rP.thresh))))
            numPerm = sum(~isnan(svm.rP.hitRate),1);
            nonEmptyPermInd = find(numPerm);
            xPos = svm.pP.featLevList(nonEmptyPermInd);
            yPos = ones(size(xPos))*0.05;
            numPerm = numPerm(nonEmptyPermInd);
            for i = 1:length(numPerm)
                text(xPos(i),yPos(i),num2str(numPerm(i)))
            end
        end
    end
    
    
    
    %% Plot binomial distribution
%     if ~isfield(svm.p,'doSVM') || svm.p.doSVM
%         for i = 1:length(svm.p.binom.hitRate)
%             if isfield(svm.r,'hits')
%                 plot(svm.p.featLevList,svm.p.binom.hitRate(i)*svm.p.nObs*ones(size(svm.p.featLevList)),':r')
%             else
                alpha = 0.05;
                negThresh = binoinv(alpha,svm.p.nObs,0.5)/svm.p.nObs;
                posThresh = binoinv(1-alpha,svm.p.nObs,0.5)/svm.p.nObs;
                plot(svm.p.featLevList,posThresh*ones(size(svm.p.featLevList)),':r')
                plot(svm.p.featLevList,negThresh*ones(size(svm.p.featLevList)),':r')
                alpha = 0.025;
                negThresh = binoinv(alpha,svm.p.nObs,0.5)/svm.p.nObs;
                posThresh = binoinv(1-alpha,svm.p.nObs,0.5)/svm.p.nObs;
                plot(svm.p.featLevList,posThresh*ones(size(svm.p.featLevList)),':r')
                plot(svm.p.featLevList,negThresh*ones(size(svm.p.featLevList)),':r')
%             end
%         end
%     end
    
    %% Plot actual classifaction accuracies
    if exist('testParam','var') && ~isempty(testParam)
        svm.r.hitRate = svm.(testParam).hitRate(:,:,testParamInd);
        svm.p.(testParam) = svm.(testParam).(testParam)(testParamInd);
    end
    
%     if ~isfield(svm.p,'doSVM') || svm.p.doSVM
        if size(svm.r.hitRate,1)>1
            if isfield(svm,'pP') && isfield(svm.pP,'curAltFeatLev')
                curInd = find(svm.pP.curAltFeatLev==svm.p.featLevList)-3:find(svm.pP.curAltFeatLev==svm.p.featLevList)+3;
                if curInd(1)<=0
                    curInd = curInd + abs(curInd(1)) + 1;
                end
                if any(curInd>length(svm.p.featLevList))
                    curInd(curInd>length(svm.p.featLevList)) = [];
                end
                errorbar(svm.p.featLevList(curInd),mean(svm.r.hitRate(:,curInd),1),std(svm.r.hitRate(:,curInd),[],1)./sqrt(size(svm.r.hitRate,1)),'b')
            else
                switch svm.p.algorithm
                    case 'runSVM'
                        errorbar(svm.p.featLevList,mean(svm.r.hitRate,1),std(svm.r.hitRate,[],1)./sqrt(size(svm.r.hitRate,1)),'b')
                    case 'runSVM_RFE'
                        errorbar(fliplr(svm.p.featLevList),mean(svm.r.hitRate,1),std(svm.r.hitRate,[],1)./sqrt(size(svm.r.hitRate,1)),'b')
                        set(gca, 'xdir','reverse')
                end
            end
        else
            if size(svm.r.hitRate,2)==1
                plot(svm.p.featLevList,svm.r.hitRate,'*b')
            else
                plot(svm.p.featLevList,svm.r.hitRate,'b')
            end
        end
%     elseif svm.p.doCrossVal
%         errorbar(svm.p.featLevList,mean(svm.xValCorr.Rsquared,1),std(svm.xValCorr.Rsquared,[],1))
%         ylim([mean(mean(svm.xValCorr.Rsquared(svm.xValCorr.Rsquared>-inf),1)) 0])
%     end
    
    %% Label and limit
    if isfield(svm.p,'dataFileOut')
        titleStr = svm.p.dataFileOut;
    else
        if strcmp(svm.p.kernel,'polynomial')
            titleStr = [svm.p.subj '; ' svm.p.curParam '; ' svm.p.kernel '_' num2str(svm.p.polyorder)];
        else
            titleStr = [svm.p.subj '; ' svm.p.curParam '; ' svm.p.kernel];
        end
        if isfield(svm.p,'C')
            titleStr = [titleStr '; C=' num2str(svm.p.C)];
        end
        if isfield(svm.p,'k')
            titleStr = [titleStr '; k=' num2str(svm.p.k)];
        end
        if isfield(svm.p,'crossValType')
            titleStr = [titleStr '; ' svm.p.crossValType];
        end
        if isfield(svm.p,'repeat')
            titleStr = [titleStr '; rep=' num2str(svm.p.repeat)];
        end
        if isfield(svm.p,'dataType')
            titleStr = [titleStr '; ' svm.p.dataType];
        end
        if isfield(svm.p,'fitType')
            titleStr = [titleStr '; ' svm.p.fitType];
        end
        if isfield(svm.p,'regressSession') && svm.p.regressSession
            titleStr = [titleStr '; sessionRegOut'];
        end
        if isfield(svm.p,'pca') && svm.p.pca
            titleStr = [titleStr '; onPCA' num2str(svm.p.pca)];
        end
        if isfield(svm.p,'docICA') && svm.p.docICA
            titleStr = [titleStr '; featTrans_' char(svm.p.curcICAalgo)];
        end
        if isfield(svm.p,'rotateData') && svm.p.rotateData
            titleStr = [titleStr '; rot' num2str(svm.p.rotateData/pi*180)];
        end
        if isfield(svm.p,'algorithm')
            titleStr = [titleStr '; ' char(svm.p.algorithm)];
        end
    end
    
    
    if includeNonParamDist
        if isfield(svm.rP,'hits')
            titleStr = [titleStr '; ' num2str(size(svm.rP.hits,1)) 'perm'];
        else
            titleStr = [titleStr '; ' num2str(size(svm.rP.hitRate,1)) 'perm'];
        end
    end
    
    all_ = strfind(titleStr,'_');
    for i = length(all_):-1:1
        titleStr = [titleStr(1:all_(i)) titleStr(all_(i):end)];
    end
    titleStr_tmp = {titleStr(1:end/2); titleStr(end/2+1:end)};
    title(titleStr_tmp)
    
    xlabel('nVoxel included (sorted)')
    if oldVersion
        ylim([0 svm.p.nObs]);
    else
        ylim([0 1]);
    end
    
    if isfield(svm,'pP') && isfield(svm.pP,'curAltFeatLev')
        xlims = sort([svm.p.featLevList(curInd(1)) svm.p.featLevList(curInd(end))]);
        xlims(1) = xlims(1)-1;
        xlims(2) = xlims(2)+1;
        xlim(xlims)
    end
    
end






%% Plot Xval correlation
if isfield(svm.p,'doCrossVal') && svm.p.doCrossVal
    %% Initiate plot
    if ~exist('h2','var') || isempty(h2)
        if plotIt
            h2 = figure('windowStyle','docked');
        else
            h2 = figure('windowStyle','docked','visible','off');
        end
    else
        if plotIt
            figure(h2); hold off;
        else
            figure(h2,'visible','off'); hold off;
        end
        
    end
    
    oldVersion = 0;
    if ~isfield(svm.p,'featLevList')
        svm.p.featLevList = 1:size(svm.r.hits,2);
    end
%     if ~isfield(svm.p.binom,'hitRate')
%         oldVersion = 1;
%         svm.p.binom.hitRate = svm.p.binom.hits;
%     end
%     if ~isfield(svm.r,'hitRate')
%         svm.r.hitRate = svm.r.hits;
%     end
    
    includeNonParamDist = 0;
    if isfield(svm,'xValCorrP')
        includeNonParamDist = 1;
    end
    
    
%     if oldVersion
%         plot(svm.p.featLevList,svm.p.nObs/2*ones(size(svm.p.featLevList)),'-k'); hold on;
%     else
%         plot(svm.p.featLevList,0.5*ones(size(svm.p.featLevList)),'-k'); hold on;
%     end
    
    %% Plot non-parametric distribution
    if includeNonParamDist
        if ~isfield(svm.pP,'featLevList')
            svm.pP.featLevList = 1:size(svm.rP.hits,2);
        end
%         if ~isfield(svm.pP.binom,'hitRate')
%             oldVersion = 1;
%             svm.pP.binom.hitRate = svm.pP.binom.hits;
%         end
%         if ~isfield(svm.rP,'hitRate')
%             svm.rP.hitRate = svm.rP.hits;
%         end
        tmph = plot(0,1); delete(tmph); hold on
        for i = 1:size(svm.xValCorrP.thresh,1)
            plot(svm.pP.featLevList,svm.xValCorrP.thresh(i,:),'--r')
        end
        plot(svm.pP.featLevList,mean(svm.xValCorrP.Rsquared,1),'-r')
        
        if ((size(svm.xValCorrP.thresh,2)==1) || (size(svm.xValCorrP.thresh,2)~=size(svm.xValCorr.Rsquared,2)))
            numPerm = sum(~isnan(svm.xValCorrP.Rsquared),1);
            nonEmptyPermInd = find(numPerm);
            xPos = svm.pP.featLevList(nonEmptyPermInd);
            yPos = ones(size(xPos))*0.05;
            numPerm = numPerm(nonEmptyPermInd);
            for i = 1:length(numPerm)
                text(xPos(i),yPos(i),num2str(numPerm(i)))
            end
        end
    end
    
    
    
%     %% Plot binomial distribution
%     if ~isfield(svm.p,'doSVM') || svm.p.doSVM
%         for i = 1:length(svm.p.binom.hitRate)
%             if isfield(svm.r,'hits')
%                 plot(svm.p.featLevList,svm.p.binom.hitRate(i)*svm.p.nObs*ones(size(svm.p.featLevList)),':r')
%             else
%                 plot(svm.p.featLevList,svm.p.binom.hitRate(i)*ones(size(svm.p.featLevList)),':r')
%             end
%         end
%     end
    
    %% Plot actual classifaction accuracies
%     if exist('testParam','var')
%         svm.r.hitRate = svm.(testParam).hitRate(:,:,testParamInd);
%         svm.p.(testParam) = svm.(testParam).(testParam)(testParamInd);
%     end
    
%     if ~isfield(svm.p,'doSVM') || svm.p.doSVM
%         if size(svm.r.hitRate,1)>1
%             if isfield(svm,'pP') && isfield(svm.pP,'curAltFeatLev')
%                 curInd = find(svm.pP.curAltFeatLev==svm.p.featLevList)-3:find(svm.pP.curAltFeatLev==svm.p.featLevList)+3;
%                 if curInd(1)<=0
%                     curInd = curInd + abs(curInd(1)) + 1;
%                 end
%                 if any(curInd>length(svm.p.featLevList))
%                     curInd(curInd>length(svm.p.featLevList)) = [];
%                 end
%                 errorbar(svm.p.featLevList(curInd),mean(svm.r.hitRate(:,curInd),1),std(svm.r.hitRate(:,curInd),[],1),'b')
%             else
%                 errorbar(svm.p.featLevList,mean(svm.r.hitRate,1),std(svm.r.hitRate,[],1),'b')
%             end
%         else
%             if size(svm.r.hitRate,2)==1
%                 plot(svm.p.featLevList,svm.r.hitRate,'*b')
%             else
%                 plot(svm.p.featLevList,svm.r.hitRate,'b')
%             end
%         end
%     elseif svm.p.doCrossVal
        errorbar(svm.p.featLevList,mean(svm.xValCorr.Rsquared,1),std(svm.xValCorr.Rsquared,[],1)./sqrt(size(svm.xValCorr.Rsquared,1)))
%         ylim([mean(mean(svm.xValCorr.Rsquared(svm.xValCorr.Rsquared>-inf),1)) 0])
%     end
    
    %% Label and limit
%     if isfield(svm.p,'dataFileOut')
        titleStr = svm.p.dataFileOut;
%     else
%         if strcmp(svm.p.kernel,'polynomial')
%             titleStr = [svm.p.subj '; ' svm.p.curParam '; ' svm.p.kernel '_' num2str(svm.p.polyorder)];
%         else
%             titleStr = [svm.p.subj '; ' svm.p.curParam '; ' svm.p.kernel];
%         end
%         if isfield(svm.p,'C')
%             titleStr = [titleStr '; C=' num2str(svm.p.C)];
%         end
%         if isfield(svm.p,'k')
%             titleStr = [titleStr '; k=' num2str(svm.p.k)];
%         end
%         if isfield(svm.p,'crossValType')
%             titleStr = [titleStr '; ' svm.p.crossValType];
%         end
%         if isfield(svm.p,'repeat')
%             titleStr = [titleStr '; rep=' num2str(svm.p.repeat)];
%         end
%         if isfield(svm.p,'dataType')
%             titleStr = [titleStr '; ' svm.p.dataType];
%         end
%         if isfield(svm.p,'fitType')
%             titleStr = [titleStr '; ' svm.p.fitType];
%         end
%         if isfield(svm.p,'regressSession') && svm.p.regressSession
%             titleStr = [titleStr '; sessionRegOut'];
%         end
%         if isfield(svm.p,'pca') && svm.p.pca
%             titleStr = [titleStr '; onPCA' num2str(svm.p.pca)];
%         end
%         if isfield(svm.p,'docICA') && svm.p.docICA
%             titleStr = [titleStr '; featTrans_' char(svm.p.curcICAalgo)];
%         end
%         if isfield(svm.p,'rotateData') && svm.p.rotateData
%             titleStr = [titleStr '; rot' num2str(svm.p.rotateData/pi*180)];
%         end
%         if isfield(svm.p,'algorithm')
%             titleStr = [titleStr '; ' char(svm.p.algorithm)];
%         end
%     end
    
    
    if includeNonParamDist
        titleStr = [titleStr '; ' num2str(size(svm.xValCorrP.Rsquared,1)) 'perm'];
    end
    
    all_ = strfind(titleStr,'_');
    for i = length(all_):-1:1
        titleStr = [titleStr(1:all_(i)) titleStr(all_(i):end)];
    end
    titleStr_tmp = {titleStr(1:end/2); titleStr(end/2+1:end)};
    title(titleStr_tmp)
    
    xlabel('nVoxel included (sorted)')
    ylabel('Cross-validated R2')
%     if oldVersion
%         ylim([0 svm.p.nObs]);
%     else
%         ylim([0 1]);
%     end
    
    if isfield(svm,'pP') && isfield(svm.pP,'curAltFeatLev')
        xlims = sort([svm.p.featLevList(curInd(1)) svm.p.featLevList(curInd(end))]);
        xlims(1) = xlims(1)-1;
        xlims(2) = xlims(2)+1;
        xlim(xlims)
    end
    
end
if ~exist('h2','var')
    h2 = [];
end
if ~exist('h','var')
    h = [];
end