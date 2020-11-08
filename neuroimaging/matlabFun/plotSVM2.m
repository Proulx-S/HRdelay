function [h,h2] = plotSVM2(svm,h,testParam,testParamInd,h2,nonParamDist,plotIt)

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
    
    
    switch nonParamDist
        case 'actualMean'
            tmp = squeeze(nanmean(svm.r.hitRate,1))';
        case 'negThresh'
            tmp = squeeze(nanmean(svm.rP.thresh(2,:,:),1))';
        case 'posThresh'
            tmp = squeeze(nanmean(svm.rP.thresh(3,:,:),1))';
        case 'actualThresh'
            tmp = squeeze(nanmean(svm.r.hitRate,1))';
            tmpNeg = squeeze(nanmean(svm.rP.thresh(2,:,:),1))';
            tmpPos = squeeze(nanmean(svm.rP.thresh(3,:,:),1))';
            tmp(tmp>tmpNeg & tmp<tmpPos) = nan;
        case 'permMean'
            tmp = squeeze(nanmean(svm.rP.hitRate,1))';
        case 'actualBinomDistThresh'
            tmp = binomialSigThresh2(svm)';
%             tmp = squeeze(nanmean(svm.r.hitRate,1))';
%             tmpNeg = binoinv(0.05,svm.p.nObs,0.5)/svm.p.nObs;
%             tmpPos = binoinv(0.95,svm.p.nObs,0.5)/svm.p.nObs;
%             tmp(tmp>tmpNeg & tmp<tmpPos) = nan;    
    end
    dim2Name = 'fun';
    dim2Lev = svm.p.featLevList;
    dim1Name = 'feat';
    dim1Lev = svm.p.featLevList;
    
    tmp = flipud(tmp);
    dim1Lev = fliplr(dim1Lev);
    
    imagesc(tmp,[0 1]); colorbar
    xlabel(dim2Name)
    ylabel(dim1Name)
    
    
    XTick = [1 5:5:size(svm.r.hitRate,3)]; if XTick(end)~=size(svm.r.hitRate,3); XTick(end+1) = size(svm.r.hitRate,3); end;
%     hText = xticklabel_rotate(XTick,90,strread(num2str(svm.p.featLevList(XTick)),'%s'));
    set(gca,'XTick',XTick)
    if isfield(svm.p,'funLevList')
        set(gca,'XTickLabel',svm.p.funLevList(XTick))
    else
        set(gca,'XTickLabel',svm.p.featLevList(XTick))
    end
    YTick = 1:size(svm.r.hitRate,2);
%     YTick = [1 5:5:size(svm.r.hitRate,2)]; if YTick(end)~=size(svm.r.hitRate,2); YTick(end+1) = size(svm.r.hitRate,2); end;
    set(gca,'YTick',YTick)
    set(gca,'YTickLabel',fliplr(svm.p.featLevList(YTick)))
    
    

    
    if isfield(svm.p,'dataFileOut')
        titleStr = svm.p.dataFileOut;
    else
        titleStr = [];
    end
    
    
    if ~strcmp(nonParamDist,'actualMean') && ~strcmp(nonParamDist,'actualBinomDistThresh')
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