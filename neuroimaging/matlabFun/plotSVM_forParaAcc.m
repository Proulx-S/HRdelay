function h = plotSVM_forParaAcc(svm,h,plotIt,absFlag)
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

if length(size(svm.rTe{1}.acc))>3
    [acc,info] = getAccGroup(svm.rTe,absFlag);
    %Collapse different sessions
    for i = 1:length(svm.pAll)
        subjList{i} = svm.pAll{i}.subj;
        sessList{i} = svm.pAll{i}.preSelection;
    end
    sess1ind = strmatch('sess1',sessList);
    sess2ind = strmatch('sess2',sessList);
    for i = 1:length(sess1ind)
        sess1cur = sess1ind(i);
        tmp = strmatch(subjList(sess1ind(i)),subjList);
        sess2cur = tmp(ismember(tmp,sess2ind));
        
        if isfield(svm.rTe{1},'n')
            n(i,1) = svm.rTe{1}.n(sess1cur);
            n(i,2) = svm.rTe{1}.n(sess2cur);
        end
        acc2(:,:,:,i,1) = acc(:,:,:,sess1cur);
        acc2(:,:,:,i,2) = acc(:,:,:,sess2cur);
    end
    if isfield(svm.rTe{1},'n')
        for i = 1:length(svm.rTe)
            svm.rTr{i}.n = sum(n,2);
            svm.rTe{i}.n = sum(n,2);
        end
    end
    acc = mean(acc2,5);
    acc = squeeze(mean(acc,3));
else
    acc = []; info = [];
    for i = 1:length(svm.rTe)
        [accTmp,infoTmp] = getAcc(svm.rTe(i),absFlag);
        acc = cat(2,accTmp,acc);
        info = cat(2,infoTmp,info);
    end
end

accMean = mean(acc,3);
if isfield(svm.rTe{1},'n')
    accStd = std(acc,[],3)./sqrt(length(svm.rTe{1}.n));
else
    accStd = std(acc,[],3);
end

if absFlag
    alpha = 0.025;
else
    alpha = 0.05;
end
if isfield(svm.rTe{1},'n')
    negThresh = binoinv(alpha,sum(svm.rTe{1}.n),0.5)/sum(svm.rTe{1}.n)*100;
    posThresh = binoinv(1-alpha,sum(svm.rTe{1}.n),0.5)/sum(svm.rTe{1}.n)*100;
else
    negThresh = binoinv(alpha,svm.p.nObs,0.5)/svm.p.nObs*100;
    posThresh = binoinv(1-alpha,svm.p.nObs,0.5)/svm.p.nObs*100;
end


if plotPerm
    accPerm = [];
    for i = 1:length(svm.perm.rTe)
        [accTmp,~] = getAcc(svm.perm.rTe(i),absFlag);
        accPerm = cat(2,accTmp,accPerm);
    end
    accPermMean = mean(accPerm,3);
    accPerm5perc = prctile(accPerm,[alpha*100 (1-alpha)*100],3);
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

if absFlag
    tmpStr = 'absolute accuracy';
else
    tmpStr = 'accuracy';
end
if isfield(svm.rTe{1},'n')
    ylabel([tmpStr ' (I=sem)'])
else
    ylabel([tmpStr ' (I=std)'])
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
if ~absFlag
    plot(Xlim,[negThresh negThresh],'r')
end

if plotPerm
    if ~iscell(xOffset) && xOffset==0
        errorbar(Xtick,accPermMean,accPermMean-accPerm5perc(:,:,1),accPerm5perc(:,:,2)-accPermMean,'.r')
    else
        for i = 1:length(xOffset)
            errorbar(Xtick+xOffset{i},accPermMean(:,i),accPermMean(:,i)-accPerm5perc(:,i,1),accPerm5perc(:,i,2)-accPermMean(:,i),'.r')
        end
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

function [acc,info] = getAcc(r,absFlag)
sz = size(r{1}.acc);
if length(sz)<3
    sz(3) = 0;
end
for i = 1:length(r)
    subLength(i) = length(r{i}.subLevel);
end
acc = nan(length(r),max(subLength)+1,size(r{1}.acc,length(sz)));
% acc = nan(length(r),max(subLength)+1,size(r{1}.acc));
info = cell(length(r),max(subLength)+1);
for i = 1:length(r)
    acc(i,1,:) = recomputeAcc(r{i},absFlag);
    if isfield(r{i},'info')
        info{i,1} = r{i}.info;
    else
        info{i,1} = '';
    end
    for ii = 1:length(r{i}.subLevel)
        acc(i,1+ii,:) = recomputeAcc(r{i}.subLevel{ii},absFlag);
        if isfield(r{i}.subLevel{ii},'info')
            info{i,1+ii} = r{i}.subLevel{ii}.info;
        else
            info{i,1+ii} = '';
        end
    end
end

function [acc,info] = getAccGroup(r,absFlag)
sz = size(r{1}.acc);
acc = nan(length(r),length(r{1}.subLevel)+1,size(r{1}.acc,3),size(r{1}.acc,4));
info = cell(length(r),length(r{1}.subLevel)+1);
for i = 1:length(r)
    acc(i,1,:,:) = recomputeAccGroup(r{i},absFlag);
    if isfield(r{i},'info')
        info{i,1} = r{i}.info;
    else
        info{i,1} = '';
    end
    for ii = 1:length(r{i}.subLevel)
        acc(i,1+ii,:,:) = recomputeAccGroup(r{i}.subLevel{ii},absFlag);
        if isfield(r{1}.subLevel{ii},'info')
            info{i,1+ii} = r{1}.subLevel{ii}.info;
        else
            info{i,1+ii} = '';
        end
    end
end
