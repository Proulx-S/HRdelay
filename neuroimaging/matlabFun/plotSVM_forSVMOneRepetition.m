function h = plotSVM_forSVMOneRepetition(svm,h,plotIt)
if isfield(svm,'rP')
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



%% Compute data to plot
allFields = fields(svm.r.newRes);
accInd = sort([find(strcmp(allFields,'hitRate')); find(~cellfun('isempty',regexp(allFields,'hitRateL.')))]);
infoInd = sort([find(strcmp(allFields,'info')); find(~cellfun('isempty',regexp(allFields,'infoL.')))]);
acc = nan(size(svm.r.newRes.info,2),length(accInd),size(svm.r.newRes.(allFields{accInd(1)}),1)); %analysis x subLevel x rep
if plotPerm
    accP = nan(size(svm.rP.newRes.info,2),length(accInd),size(svm.rP.newRes.(allFields{accInd(1)}),1)); %analysis x subLevel x rep
end
info = cell(size(svm.r.newRes.info,2),length(accInd));
info2 = cell(size(svm.r.newRes.info,2),length(accInd));

for ii = 1:length(accInd)
    acc(:,ii,:) = permute(svm.r.newRes.(allFields{accInd(ii)}),[2 3 1]);
    if plotPerm
        accP(:,ii,:) = permute(svm.rP.newRes.(allFields{accInd(ii)}),[2 3 1]);
    end
    info(:,ii) = svm.r.newRes.(allFields{infoInd(ii)})(1,:)';
    info2(:,ii) = strcat(svm.r.newRes.(allFields{infoInd(1)})(1,:)',repmat({'__'},size(acc,1),1),svm.r.newRes.(allFields{infoInd(ii)})(1,:)');
end
% % reshape(info',1,numel(info))
% ~cellfun('isempty',strfind(info,'n/a'))

accMean = nanmean(acc,3);
accStd = nanstd(acc,[],3);

negThresh = binoinv(0.05,svm.p.nObs,0.5)/svm.p.nObs;
posThresh = binoinv(1-0.05,svm.p.nObs,0.5)/svm.p.nObs;

if plotPerm
    %% Compute non-param thresholds
    negThresh_nonParam = prctile(accP,5,3);
    posThresh_nonParam = prctile(accP,95,3);
    mean_nonParam = nanmean(accP,3);
end

%% Plot it
XTickLabel = info(:,1);
Xtick = 1:length(XTickLabel);
hb = bar(accMean);
set(gca,'XTickLabelMode','manual');
set(gca,'Xtick',Xtick);
set(gca,'XTickLabel',XTickLabel);
set(gca,'XTickLabelRotation',45)

ylabel('accuracy (I=std)')
ylim([0 1])
hold on
xOffset = get(hb,'xOffset');
for i = 1:length(xOffset)
    errorbar(Xtick+xOffset{i},accMean(:,i),accStd(:,i),'.k')
end

Xlim = get(gca,'Xlim');
plot(Xlim,[0.5 0.5],'k')
plot(Xlim,[posThresh posThresh],'r')


if plotPerm
    for i = 1:length(xOffset)
        errorbar(Xtick+xOffset{i},mean_nonParam(:,i),negThresh_nonParam(:,i)-mean_nonParam(:,i),posThresh_nonParam(:,i)-mean_nonParam(:,i),'.r')
    end
end


if isfield(svm.p,'dataFileOut')
    titleStr = svm.p.dataFileOut;
else
    titleStr = [];
end

all_ = strfind(titleStr,'_');
for i = length(all_):-1:1
    titleStr = [titleStr(1:all_(i)) titleStr(all_(i):end)];
end
tmpLength = length(titleStr);
titleStr_tmp = {titleStr(1:round(tmpLength/2)); titleStr(round(tmpLength/2)+1:end)};
if plotPerm
    if strcmp(svm.pP.algorithm,'runSVMOneRepetition')
        titleStr_tmp = [titleStr_tmp; {[num2str(size(svm.rP.newRes.hitRate,1)) 'perm']}];
    else
        titleStr_tmp = [titleStr_tmp; {[num2str(size(svm.rP.hitRate,1)) 'perm']}];
    end
end
title(titleStr_tmp,'interpret','none')

legend([{'all'} info(end,2:end)],'location','best','interpreter','none')



