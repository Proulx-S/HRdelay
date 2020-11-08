function h = plotSVM_forParaDistT(svm,h,plotIt,absFlag)
if isfield(svm,'perm')
    plotPerm = 1;
else
    plotPerm = 0;
end
if ~exist('plotIt','var') || isempty(plotIt)
    plotIt = 1;
end
if ~exist('absFlag','var') || isempty(absFlag)
    absFlag = 0;
end


%% Plot SVM
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

[dist,info] = getDist(svm.rTe,absFlag);
distMean = mean(dist,3);
if isfield(svm.rTe{1},'n')
    distStd = std(dist,[],3)./sqrt(length(svm.rTe{1}.n));
else
    distStd = std(dist,[],3);
end

% if isfield(svm.rTe{1},'n')
%     negThresh = binoinv(0.05,sum(svm.rTe{1}.n),0.5)/sum(svm.rTe{1}.n)*100;
%     posThresh = binoinv(1-0.05,sum(svm.rTe{1}.n),0.5)/sum(svm.rTe{1}.n)*100;
% else
%     negThresh = binoinv(0.05,svm.p.nObs,0.5)/svm.p.nObs*100;
%     posThresh = binoinv(1-0.05,svm.p.nObs,0.5)/svm.p.nObs*100;
% end


if plotPerm
    [distPerm,info] = getDist(svm.perm.rTe,absFlag);
    distPermMean = mean(distPerm,3);
    distPerm5perc = prctile(distPerm,[5 95],3);
end





%% Plot it
if size(info,1)==1
    XTickLabel = info;
else
    XTickLabel = info(:,1);
end
for i = 1:length(XTickLabel)
    XTickLabel{i} = strjoin(strsplit(XTickLabel{i},'_'),'__');
end
Xtick = 1:length(XTickLabel);
hb = bar(distMean);
set(gca,'XTickLabelMode','manual');
set(gca,'Xtick',Xtick);
set(gca,'XTickLabel',XTickLabel);
set(gca,'XTickLabelRotation',45)

if isfield(svm.rTe{1},'n')
    ylabel('accuracy (I=sem)')
else
    ylabel('accuracy (I=std)')
end
% ylim([0 100])
hold on
xOffset = get(hb,'xOffset');
if ~iscell(xOffset) && xOffset==0
    errorbar(Xtick,distMean,distStd,'.k')
else
    for i = 1:length(xOffset)
        errorbar(Xtick+xOffset{i},distMean(:,i),distStd(:,i),'.k')
    end
end

% Xlim = get(gca,'Xlim');
% plot(Xlim,[50 50],'k')
% plot(Xlim,[posThresh posThresh],'r')
% plot(Xlim,[negThresh negThresh],'r')


if plotPerm
    if ~iscell(xOffset) && xOffset==0
        errorbar(Xtick,distMean,distStd,'.k')
        errorbar(Xtick,distPermMean(:,i),distPermMean(:,i)-distPerm5perc(:,i,1),distPerm5perc(:,i,2)-distPermMean(:,i),'.r')
    else
        for i = 1:length(xOffset)
            errorbar(Xtick+xOffset{i},distPermMean(:,i),distPermMean(:,i)-distPerm5perc(:,i,1),distPerm5perc(:,i,2)-distPermMean(:,i),'.r')
        end
    end

    
    
    for i = 1:length(xOffset)
        errorbar(Xtick+xOffset{i},distPermMean(:,i),distPermMean(:,i)-distPerm5perc(:,i,1),distPerm5perc(:,i,2)-distPermMean(:,i),'.r')
    end
end


if isfield(svm.p,'dataFileOut')
    titleStr = svm.p.dataFileOut;
else
    titleStr = [];
end

% all_ = strfind(titleStr,'_');
% for i = length(all_):-1:1
%     titleStr = [titleStr(1:all_(i)) titleStr(all_(i):end)];
% end
tmpLength = length(titleStr);
titleStr_tmp = {titleStr(1:round(tmpLength/2)); titleStr(round(tmpLength/2)+1:end)};
if plotPerm
    titleStr_tmp = [titleStr_tmp; {[num2str(size(distPerm,3)) 'perm']}];
end
title(titleStr_tmp,'interpret','none')

if size(distMean,1)>1
    legend(info(1,:)','location','bestOutside','interpreter','none')
end

function [dist,info] = getDist(r,absFlag)
sz = size(r{1}.acc);
dist = nan(length(r),length(r{1}.subLevel)+1,size(r{1}.acc,length(sz)));
info = cell(length(r),length(r{1}.subLevel)+1);
for i = 1:length(r)
    dist(i,1,:) = computeDist(r{i},absFlag);
    if isfield(r{i},'info')
        info{i,1} = r{i}.info;
    else
        info{i,1} = '';
    end
    for ii = 1:length(r{i}.subLevel)
        dist(i,1+ii,:) = computeDist(r{i}.subLevel{ii},absFlag);
        if isfield(r{i},'info')
        info{i,1+ii} = r{1}.subLevel{ii}.info;
        else
            info{i,1+ii} = '';
        end
    end
end