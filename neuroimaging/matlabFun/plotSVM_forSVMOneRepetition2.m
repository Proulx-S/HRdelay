function h = plotSVM_forSVMOneRepetition2(svm,h,plotIt)
if isfield(svm,'perm')
    plotPerm = 1;
else
    plotPerm = 0;
end

% if isfield(svm,'plotIt') && isnumeric(svm.plotIt) && svm.plotIt
%     plotIt = svm.plotIt;
% else
%     plotIt = 1;
% end
if ~exist('plotIt','var')
    plotIt = 1;
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
acc = nan(length(svm.rTe),length(svm.rTe{1}.subLevel)+1,length(svm.rTe{1}.acc));
info = cell(length(svm.rTe),length(svm.rTe{1}.subLevel)+1);
for i = 1:length(svm.rTe)
    acc(i,1,:) = svm.rTe{i}.acc;
    info{i,1} = svm.rTe{i}.info;
    for ii = 1:length(svm.rTe{i}.subLevel)
        acc(i,1+ii,:) = svm.rTe{i}.subLevel{ii}.acc;
        info{i,1+ii} = svm.rTe{1}.subLevel{ii}.info;
    end
end
accMean = mean(acc,3);
if isfield(svm.rTe{1},'n')
    accStd = std(acc,[],3)./sqrt(length(svm.rTe{1}.n));
else
    accStd = std(acc,[],3);
end

if isfield(svm.rTe{1},'n')
    negThresh = binoinv(0.05,sum(svm.rTe{1}.n),0.5)/sum(svm.rTe{1}.n)*100;
    posThresh = binoinv(1-0.05,sum(svm.rTe{1}.n),0.5)/sum(svm.rTe{1}.n)*100;
else
    negThresh = binoinv(0.05,svm.p.nObs,0.5)/svm.p.nObs*100;
    posThresh = binoinv(1-0.05,svm.p.nObs,0.5)/svm.p.nObs*100;
end


if plotPerm
    accPerm = nan(length(svm.perm.rTe),length(svm.perm.rTe{1}.subLevel)+1,length(svm.perm.rTe{1}.acc));
    for i = 1:length(svm.perm.rTe)
        accPerm(i,1,:) = svm.perm.rTe{i}.acc;
        for ii = 1:length(svm.perm.rTe{i}.subLevel)
            accPerm(i,1+ii,:) = svm.perm.rTe{i}.subLevel{ii}.acc;
        end
    end
    accPermMean = mean(accPerm,3);
    accPerm5perc = prctile(accPerm,[5 95],3);
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
hb = bar(accMean);
set(gca,'XTickLabelMode','manual');
set(gca,'Xtick',Xtick);
set(gca,'XTickLabel',XTickLabel);
set(gca,'XTickLabelRotation',45)

if isfield(svm.rTe{1},'n')
    ylabel('accuracy (I=sem)')
else
    ylabel('accuracy (I=std)')
end
ylim([0 100])
hold on
xOffset = get(hb,'xOffset');
if ~iscell(xOffset) && xOffset==0
    errorbar(Xtick,accMean,accStd,'.k')
else
    for i = 1:length(xOffset)
        errorbar(Xtick+xOffset{i},accMean(:,i),accStd(:,i),'.k')
    end
end

Xlim = get(gca,'Xlim');
plot(Xlim,[50 50],'k')
plot(Xlim,[posThresh posThresh],'r')
plot(Xlim,[negThresh negThresh],'r')


if plotPerm
    for i = 1:length(xOffset)
        errorbar(Xtick+xOffset{i},accPermMean(:,i),accPermMean(:,i)-accPerm5perc(:,i,1),accPerm5perc(:,i,2)-accPermMean(:,i),'.r')
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
    titleStr_tmp = [titleStr_tmp; {[num2str(size(accPerm,3)) 'perm']}];
end
title(titleStr_tmp,'interpret','none')

if size(accMean,1)>1
    legend(info(1,:)','location','bestOutside','interpreter','none')
end

