function h = plotSVM_forParaDist(svm,h,plotIt,absFlag,tFlag)
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
if ~exist('tFlag','var') || isempty(tFlag)
    tFlag = 0;
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
    [dist,info] = getDistGroup(svm.rTe,absFlag,tFlag);
    
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
        dist2(:,:,:,i,1) = dist(:,:,:,sess1cur);
        dist2(:,:,:,i,2) = dist(:,:,:,sess2cur);
    end
    if isfield(svm.rTe{1},'n')
        for i = 1:length(svm.rTr)
            svm.rTr{i}.n = sum(n,2);
            svm.rTe{i}.n = sum(n,2);
        end
    end
    dist = mean(dist2,5);
    dist = squeeze(mean(dist,3));
else
    dist = []; info = [];
    for i = 1:length(svm.rTe)
        [distTmp,infoTmp] = getDist(svm.rTe(i),absFlag,tFlag);
        dist = cat(2,distTmp,dist);
        info = cat(2,infoTmp,info);
    end
%     [dist,info] = getDist(svm.rTe,absFlag,tFlag);
end
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
    [distPerm,~] = getDist(svm.perm.rTe,absFlag,tFlag);
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

if absFlag; tmpStr1 = 'absolute '; else tmpStr1 = ''; end
if tFlag; tmpStr2 = 'distT'; else tmpStr2 = 'dist'; end
if isfield(svm.rTe{1},'n')
    ylabel([tmpStr1 tmpStr2 ' (I=sem)'])
%     ylabel('accuracy (I=sem)')
else
    ylabel([tmpStr1 tmpStr2 ' (I=std)'])
%     ylabel('accuracy (I=std)')
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
        errorbar(Xtick,distPermMean,distPermMean-distPerm5perc(:,:,1),distPerm5perc(:,:,2)-distPermMean,'.r')
    else
        for i = 1:length(xOffset)
            errorbar(Xtick+xOffset{i},distPermMean(:,i),distPermMean(:,i)-distPerm5perc(:,i,1),distPerm5perc(:,i,2)-distPermMean(:,i),'.r')
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
    titleStr_tmp = [titleStr_tmp; {[num2str(size(distPerm,3)) 'perm']}];
end
title(titleStr_tmp,'interpret','none')

if size(distMean,1)>1
    legend(info(1,:)','location','bestOutside','interpreter','none')
end

function [dist,info] = getDist(r,absFlag,tFlag)
sz = size(r{1}.acc);
dist = nan(length(r),length(r{1}.subLevel)+1,size(r{1}.acc,length(sz)));
info = cell(length(r),length(r{1}.subLevel)+1);
for i = 1:length(r)
    dist(i,1,:) = computeDist(r{i},absFlag,tFlag);
    if isfield(r{i},'info')
        info{i,1} = r{i}.info;
    else
        info{i,1} = '';
    end
    for ii = 1:length(r{i}.subLevel)
        dist(i,1+ii,:) = computeDist(r{i}.subLevel{ii},absFlag,tFlag);
        if isfield(r{i},'info')
        info{i,1+ii} = r{1}.subLevel{ii}.info;
        else
            info{i,1+ii} = '';
        end
    end
end

function [dist,info] = getDistGroup(r,absFlag,tFlag)
sz = size(r{1}.acc);
dist = nan(length(r),length(r{1}.subLevel)+1,size(r{1}.acc,3),size(r{1}.acc,4));
info = cell(length(r),length(r{1}.subLevel)+1);
for i = 1:length(r)
    dist(i,1,:,:) = computeDistGroup(r{i},absFlag,tFlag);
    if isfield(r{i},'info')
        info{i,1} = r{i}.info;
    else
        info{i,1} = '';
    end
    for ii = 1:length(r{i}.subLevel)
        dist(i,1+ii,:,:) = computeDistGroup(r{i}.subLevel{ii},absFlag,tFlag);
        if isfield(r{i},'info')
        info{i,1+ii} = r{1}.subLevel{ii}.info;
        else
            info{i,1+ii} = '';
        end
    end
end