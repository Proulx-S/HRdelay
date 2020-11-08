function h = plotSVM_forSVMOneRepetition3(svm,h,plotIt)
if isfield(svm,'perm') && ~isempty(svm.perm.rTe{1}.acc)
    plotPerm = 1;
else
    plotPerm = 0;
end

% if isfield(svm,'plotIt') && isnumeric(svm.plotIt) && svm.plotIt
%     plotIt = svm.plotIt;
% else
%     plotIt = 1;
% end
if ~exist('plotIt','var') || isempty(plotIt)
    plotIt = 1;
end


%% Plot SVM
%% Initiate plot
if ~exist('h','var') || isempty(h)
    for anaInd = 1:length(svm.rTr)
        if plotIt
            h(anaInd) = figure('windowStyle','docked');
        else
            h(anaInd) = figure('visible','off');
        end
        colormap jet
    end
else
    for anaInd = 1:length(svm.rTr)
        if plotIt
            figure(h(anaInd)); hold off;
        else
            figure(h(anaInd),'visible','off'); hold off;
        end
    end
end

%% Extract data
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
accStd = std(acc,[],3);

negThresh = binoinv(0.05,svm.p.nObs,0.5)/svm.p.nObs*100;
posThresh = binoinv(1-0.05,svm.p.nObs,0.5)/svm.p.nObs*100;


%% Extract permutations
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

%% Make one figure for each analysis
accMean_orig = accMean;
accStd_orig = accStd;
info_orig = info;
if plotPerm
    accPermMean_orig = accPermMean;
    accPerm5perc_orig = accPerm5perc;
end


for anaInd = 1:size(accMean_orig,1)
    figure(h(anaInd));
    
    accMean = accMean_orig(anaInd,:);
    accStd = accStd_orig(anaInd,:);
    info = info_orig(anaInd,:);
    if plotPerm
        accPermMean = accPermMean_orig(anaInd,:);
        accPerm5perc = squeeze(accPerm5perc_orig(anaInd,:,:))';
    end

    %------------------------------------------------------
    % Special case for 'ALrev'
    tmp = [];
    for i = 1:length(svm.p.infoComb)
        tmp = [tmp strcmp(svm.p.infoComb{i},'ALrev1')];
    end
    if all(tmp)
%     if sepFig
        accMean1 = accMean(1,1); accMean(1) = [];
        accMean1 = cat(1,nan,nan);
        accMean = cat(1,accMean(1:end/2),accMean(end/2+1:end));
        accMean = cat(2,accMean1,accMean); clear accMean1
        
        accStd1 = accStd(1,1); accStd(1) = [];
        accStd1 = cat(1,nan,nan);
        accStd = cat(1,accStd(1:end/2),accStd(end/2+1:end));
        accStd = cat(2,accStd1,accStd); clear accStd1
        
        tmp = info(1,1); info(1) = [];
        tmp = cat(1,{''},{''});
        info = cat(1,info(1:end/2),info(end/2+1:end));
        info = cat(2,tmp,info); clear info1
        info{1,1} = ['before ' svm.p.infoComb{anaInd}];
        info{2,1} = ['after ' svm.p.infoComb{anaInd}];
        
        if plotPerm
            tmp = accPermMean(1,1); accPermMean(1) = [];
            tmp = cat(1,nan,nan);
            accPermMean = cat(1,accPermMean(1:end/2),accPermMean(end/2+1:end));
            accPermMean = cat(2,tmp,accPermMean);
            
            accPerm5perc = permute(accPerm5perc,[3 2 1]);
            tmp = accPerm5perc(1,1,:); accPerm5perc(:,1,:) = [];
            tmp = cat(1,tmp,tmp); tmp = nan(size(tmp));
            accPerm5perc = cat(1,accPerm5perc(:,1:end/2,:),accPerm5perc(:,end/2+1:end,:));
            accPerm5perc = cat(2,tmp,accPerm5perc);
        end
    end
    %------------------------------------------------------
    
    
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
    
    ylabel('accuracy (I=std)')
    ylim([0 100])
    hold on
    xOffset = get(hb,'xOffset');
    if ~iscell(xOffset) && xOffset==0
        errorbar(Xtick,accMean,accStd,'.k')
        if plotPerm
            errorbar(Xtick,accPermMean,accPermMean-accPerm5perc(1,:),accPerm5perc(2,:)-accPermMean,'or')
        end
    else
        for i = 1:length(xOffset)
            errorbar(Xtick+xOffset{i},accMean(:,i),accStd(:,i),'.k')
            if plotPerm
                errorbar(Xtick+xOffset{i},accPermMean(:,i),accPermMean(:,i)-accPerm5perc(:,i,1),accPerm5perc(:,i,1)-accPermMean(:,i),'or')
%                 errorbar(Xtick+xOffset{i},accPermMean(:,i),accPermMean-accPerm5perc(1,:,1),accPerm5perc(2,:,1)-accPermMean(:,:,1),'or')
%                 error('too lazy to code that')
%                 errorbar(Xtick+xOffset{i},accMean(:,i),accStd(:,i),'.k')
            end
        end
    end
    
    Xlim = get(gca,'Xlim');
    plot(Xlim,[50 50],'k')
    plot(Xlim,[posThresh posThresh],'r')
    
    
%     if plotPerm
%         for i = 1:length(xOffset)
%             errorbar(Xtick+xOffset{i},mean_nonParam(:,i),negThresh_nonParam(:,i)-mean_nonParam(:,i),posThresh_nonParam(:,i)-mean_nonParam(:,i),'.r')
%         end
%     end
    
    
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
        titleStr_tmp = [titleStr_tmp; {[num2str(length(svm.perm.rTe{1}.acc)) 'perm']}];
    end
    title(titleStr_tmp,'interpret','none')
    
    if size(accMean,1)>1
        legend(info(1,:)','location','bestOutside','interpreter','none')
    end
end

