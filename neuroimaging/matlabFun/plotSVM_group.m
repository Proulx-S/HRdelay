function h = plotSVM_group(metric,metricPerm,fieldName)
if ~exist('metricPerm','var') || isempty(metricPerm)
    plotPerm = 0;
else
    plotPerm = 1;
end



%% Plot SVM
%% Initiate plot
h = figure('windowStyle','docked');

for anaInd = 1:length(metric)
    info{anaInd,1} = metric{anaInd}.info;
    data = mean(metric{anaInd}.(fieldName),2);
    dataMean(anaInd,1) = nanmean(data,1);
    dataSem(anaInd,1) = nanstd(data,[],1)/sqrt(size(data(~isnan(data),:),1));
    for i = 1:length(metric{anaInd}.subLevel)
        info{anaInd,1+i} = metric{anaInd}.subLevel{i}.info;
        data = mean(metric{anaInd}.subLevel{i}.(fieldName),2);
        dataMean(anaInd,1+i) = nanmean(data,1);
        dataSem(anaInd,1+i) = nanstd(data,[],1)/sqrt(size(data(~isnan(data),:),1));
    end
end

switch fieldName
    case {'auc','acc','dist','distT'}
        alpha = 0.025;
    case {'accAbs','distAbs','distTabs'}
        alpha = 0.05;
end

negThresh = binoinv(alpha,sum(metric{anaInd}.n(:,1)),0.5)/sum(metric{anaInd}.n(:,1))*100;
posThresh = binoinv(1-alpha,sum(metric{anaInd}.n(:,1)),0.5)/sum(metric{anaInd}.n(:,1))*100;


if plotPerm
    for anaInd = 1:length(metricPerm)
        data = metricPerm{anaInd}.(fieldName);
        dataPermMean(anaInd,1) = mean(data,2);
        data5percPerm(anaInd,1,:) = prctile(data,[5 95],2);
        for i = 1:length(metricPerm{anaInd}.subLevel)
            data = metricPerm{anaInd}.subLevel{i}.(fieldName);
            dataPermMean(anaInd,1+i) = mean(data,2);
            data5percPerm(anaInd,1+i,:) = prctile(data,[5 95],2);
        end
    end
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
hb = bar(dataMean);
set(gca,'XTickLabelMode','manual');
set(gca,'Xtick',Xtick);
set(gca,'XTickLabel',XTickLabel);
set(gca,'XTickLabelRotation',45)


ylabel([fieldName ' (I=sem)'])
switch fieldName
    case {'auc','acc'}
        ylim([0 100])
    case 'accAbs'
            ylim([50 100])
    case {'dist','distAbs','distT','distTabs'}
    otherwise
        error('XX')
end
hold on
xOffset = get(hb,'xOffset');
if ~iscell(xOffset) && xOffset==0
    errorbar(Xtick,dataMean,dataSem,'.k')
else
    for i = 1:length(xOffset)
        errorbar(Xtick+xOffset{i},dataMean(:,i),dataSem(:,i),'.k','marker','none')
    end
end

Xlim = get(gca,'Xlim');

switch fieldName
    case 'auc'
        plot(Xlim,[50 50],'k')
    case 'acc'
        plot(Xlim,[50 50],'k')
        plot(Xlim,[posThresh posThresh],'r')
        plot(Xlim,[negThresh negThresh],'r')
    case {'dist','distT'}
        plot(Xlim,[0 0],'k')
    case {'distAbs','distTabs','accAbs'}
    otherwise
        error('XX')
end

if plotPerm
    if ~iscell(xOffset) && xOffset==0
        errorbar(Xtick,dataPermMean,dataPermMean-data5percPerm(:,:,1),data5percPerm(:,:,2)-dataPermMean,'.r')
    else
        for i = 1:length(xOffset)
            errorbar(Xtick+xOffset{i},dataPermMean(:,i),dataPermMean(:,i)-data5percPerm(:,i,1),data5percPerm(:,i,2)-dataPermMean(:,i),'.r')
        end
    end
end

titleStr = 'group summary';
if plotPerm
    titleStr = [titleStr ' ' num2str(size(metricPerm{anaInd}.acc,2)) ' perm'];
end
title(titleStr,'interpret','none')

if size(dataMean,1)>1
    for ii = 1:numel(info)
        if isempty(info{ii})
            info{ii} = 'n/a';
        end
    end
    info2 = info(1,:)';
    for iii = 2:size(info,1)
        info2 = strcat(info2,repmat({';   '},[size(info,2) 1]),info(iii,:)');
    end
    info2{1} = 'comb';
    legend(info2,'location','bestOutside','interpreter','none')
end
