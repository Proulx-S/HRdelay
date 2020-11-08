function [h,h2] = plotSVM_patAct(svm,h,testParam,testParamInd,h2,nonParamDist,plotIt)

% if isfield(svm,'plotIt') && isnumeric(svm.plotIt) && svm.plotIt
%     plotIt = svm.plotIt;
% else
%     plotIt = 1;
% end
if ~exist('plotIt','var')
    plotIt = 1;
end  

if ~exist('nonParamDist','var') || isempty(nonParamDist)
    nonParamDist = 'actualMean';
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
    
    tmp = squeeze(nanmean(nanmean(svm.r.hitRate,1),4));
    negThresh = binoinv(0.05,svm.p.nObs,0.5)/svm.p.nObs;
    posThresh = binoinv(1-0.05,svm.p.nObs,0.5)/svm.p.nObs;
    tmp(tmp>negThresh & tmp<posThresh) = nan;
    
    
    imagesc(tmp,[0 1]); colorbar
    xlabel('test')
    ylabel('train')
    Tick = 1:length(tmp);
    set(gca,'YTick',Tick)
    set(gca,'YTickLabel',svm.r.crossInfo_dim1Train)
    set(gca,'XTick',Tick)
    set(gca,'XTickLabel',svm.r.crossInfo_dim2Test)
    set(gca, 'XAxisLocation', 'top')
    set(gca, 'TickDir', 'out')
    

    
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
    title(titleStr_tmp,'interpret','none')

end







if ~exist('h2','var')
    h2 = [];
end
if ~exist('h','var')
    h = [];
end