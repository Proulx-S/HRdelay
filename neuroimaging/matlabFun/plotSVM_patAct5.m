function [h,h2] = plotSVM_patAct5(svm,h,nonParamDist,plotIt)

% if isfield(svm,'plotIt') && isnumeric(svm.plotIt) && svm.plotIt
%     plotIt = svm.plotIt;
% else
%     plotIt = 1;
% end
if ~exist('plotIt','var')
    plotIt = 1;
end  

if ~exist('nonParamDist','var') || isempty(nonParamDist)
    nonParamDist = 0;
end  

%% Plot SVM
if ~isfield(svm.p,'doSVM') || svm.p.doSVM
    %% Initiate plot
    if ~exist('h','var') || isempty(h)
        if plotIt
            h = figure('windowStyle','docked');
        else
            h = figure('visible','off');
        end
        colormap jet
    else
        if plotIt
            figure(h); hold off;
        else
            figure(h,'visible','off'); hold off;
        end
        
    end
    
    %% Compute data to plot
    if isfield(svm.r,'newRes')
        acc = svm.r.newRes.hitRate;
        accMean = squeeze(nanmean(acc,1));
        accStd = squeeze(nanstd(acc,[],1));
        svm.p.trainInfo = svm.r.newRes.info(1,:);
    else
        acc = mean(svm.r.hitRate,4);
        accMean = squeeze(nanmean(acc,1));
        accStd = squeeze(nanstd(acc,[],1));
        
        %     acc = nan(size(svm.r.hitRate,2),size(svm.r.hitRate,3),size(svm.r.hitRate,1)*size(svm.r.hitRate,4));
        %     n = 0;
        %     for i = 1:size(svm.r.hitRate,1)
        %         for ii = 1:size(svm.r.hitRate,4)
        %             n = n+1;
        %             acc(:,:,n) = squeeze(svm.r.hitRate(i,:,:,ii));
        %         end
        %     end
        %     accMean = nanmean(acc,3);
        %     accStd = nanstd(acc,[],3);
        %     accSem = accStd./sqrt(size(acc,3));
        
        accMean = diag(accMean);
        accStd = diag(accStd);
        %     accSem = diag(accSem);
    end
    
%     if isfield(svm.r,'newRes')
%         svm.p.trainInfo = svm.r.newRes.info(1,:);
%         accMean = mean(svm.r.newRes.hitRate,1)';
%         accStd = std(svm.r.newRes.hitRate,[],1)';
%     end
%     switch svm.p.infoComb
%         case 'layeredSVMwCross'
%             error('double-check that')
%             svm.p.trainInfo = [svm.p.trainInfo {[svm.p.trainInfo{1} '+' svm.p.trainInfo{2}]}];
%             accMean = [accMean; mean(svm.r.combInfo.hitRate)/100];
%             accStd = [accStd; std(svm.r.combInfo.hitRate)/100];
%             %     accSem = [accMean; std(svm.r.combInfo.hitRate)/100/sqrt(length(svm.r.combInfo.hitRate))];
%         case 'layeredSVM'
%             svm.p.trainInfo = svm.r.newRes.info(1,:);
%             accMean = mean(svm.r.newRes.hitRate,1)';
%             accStd = std(svm.r.newRes.hitRate,[],1)';
%     end
    
    
    negThresh = binoinv(0.05,svm.p.nObs,0.5)/svm.p.nObs;
    posThresh = binoinv(1-0.05,svm.p.nObs,0.5)/svm.p.nObs;
    
    if nonParamDist
        %% Compute non-param thresholds
        if ~iscell(svm.p.infoComb) && strcmp(svm.p.infoComb,'layeredSVM')
            negThresh_nonParam = svm.rP.thresh(2,:);
            posThresh_nonParam = svm.rP.thresh(3,:);
            mean_nonParam = nanmean(svm.rP.hitRate,1);
        else
            negThresh_nonParam = prctile(svm.rP.newRes.hitRate,5,1);
            posThresh_nonParam = prctile(svm.rP.newRes.hitRate,95,1);
            mean_nonParam = nanmean(svm.rP.newRes.hitRate,1);
        end
        
%         if strcmp(svm.p.infoComb,'layeredSVM')
%             negThresh_nonParam = [negThresh_nonParam; svm.rP.threshComb(2)/100];
%             posThresh_nonParam = [posThresh_nonParam; svm.rP.threshComb(3)/100];
%             mean_nonParam = [mean_nonParam; nanmean(svm.rP.hitRateComb,1)/100];
%         end
    end
    
    %% Plot it
    XTickLabel = svm.p.trainInfo;
    Xtick = 1:length(svm.p.trainInfo);
    bar(accMean)
    set(gca,'XTickLabelMode','manual');
    set(gca,'Xtick',Xtick);
    set(gca,'XTickLabel',XTickLabel);
    
    ylabel('accuracy (s.e.m.)')
    ylim([0 1])
    hold on
    errorbar(Xtick,accMean,accStd,'.')
    
    Xlim = get(gca,'Xlim');
    plot(Xlim,[0.5 0.5],'k')
    plot(Xlim,[posThresh posThresh],'r')
    
    
    if nonParamDist
        errorbar(Xtick-Xtick(1)*0.1,mean_nonParam,negThresh_nonParam-mean_nonParam,posThresh_nonParam-mean_nonParam,'ro')
    end
    
    
    if isfield(svm.p,'dataFileOut')
        titleStr = svm.p.dataFileOut;
    else
        titleStr = [];
    end
    
    
%     if ~strcmp(nonParamDist,'actualMean') && ~strcmp(nonParamDist,'actualBinomDistThresh')
%         if isfield(svm.rP,'hits')
%             titleStr = [titleStr '; ' num2str(size(svm.rP.hits,1)) 'perm'];
%         else
%             titleStr = [titleStr '; ' num2str(size(svm.rP.hitRate,1)) 'perm'];
%         end
%     end
    
    all_ = strfind(titleStr,'_');
    for i = length(all_):-1:1
        titleStr = [titleStr(1:all_(i)) titleStr(all_(i):end)];
    end
    tmpLength = length(titleStr);
    titleStr_tmp = {titleStr(1:round(tmpLength/2)); titleStr(round(tmpLength/2)+1:end)};
    if nonParamDist
        if strcmp(svm.pP.algorithm,'runSVMOneRepetition')
            titleStr_tmp = [titleStr_tmp; {[num2str(size(svm.rP.newRes.hitRate,1)) 'perm']}];
        else
            titleStr_tmp = [titleStr_tmp; {[num2str(size(svm.rP.hitRate,1)) 'perm']}];
        end
    end
    title(titleStr_tmp,'interpret','none')

end







if ~exist('h2','var')
    h2 = [];
end
if ~exist('h','var')
    h = [];
end