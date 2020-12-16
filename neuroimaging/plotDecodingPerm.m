function plotDecodingPerm(res,saveFig)
if ~exist('saveFig','var') || isempty(saveFig)
    saveFig = 0;
end

plotUpperErrorBar = 0;
lw = 1;
% yLim = [12 75]./100;
yLim = 'auto';

spaceList = fields(res)';
subjList = res.(spaceList{1}).subjList;

for spaceInd = 1:length(spaceList)
    res.(spaceList{spaceInd}).perm.acc = permute(res.(spaceList{spaceInd}).perm.acc,[2 3 1]);
    
    acc = res.(spaceList{spaceInd}).perm.acc;
    figure('WindowStyle','docked');
    yLim = [];
    for subjInd = 1:size(acc,1)
        for sessInd = 1:size(acc,2)
            subplot(size(acc,1),size(acc,2),sessInd + (subjInd-1)*2)
            histogram(acc(subjInd,sessInd,:).*100,'normalization','probability')
            yLim = [yLim; ylim];
        end
    end
    yLim = [min(yLim(:,1)) max(yLim(:,2))];
    for subjInd = 1:size(acc,1)
        for sessInd = 1:size(acc,2)
            subplot(size(acc,1),size(acc,2),sessInd + (subjInd-1)*2)
            ylim(yLim); xlim([0 100]);
            switch subjInd
                case 1
                    title(['sess' num2str(sessInd)])
                case size(acc,1)
                    xlabel('Accuracy (%correct)')
            end
            if sessInd==1
                ylabel([res.(spaceList{spaceInd}).subjList(subjInd) {'prob'}])
            end
        end
    end
    
    figure('WindowStyle','docked');
    nObs = repmat(res.(spaceList{spaceInd}).nObs,[1 1 size(acc,3)]);
    hit = round(acc.*nObs);
    size(hit)
    
    res.(spaceList{spaceInd}).perm
    histogram(acc(:).*100,'normalization','probability')
            yLim = [yLim; ylim];
        end
    end
    yLim = [min(yLim(:,1)) max(yLim(:,2))];
    for subjInd = 1:size(acc,1)
        for sessInd = 1:size(acc,2)
            subplot(size(acc,1),size(acc,2),sessInd + (subjInd-1)*2)
            ylim(yLim); xlim([0 100]);
            switch subjInd
                case 1
                    title(['sess' num2str(sessInd)])
                case size(acc,1)
                    xlabel('Accuracy (%correct)')
            end
            if sessInd==1
                ylabel([res.(spaceList{spaceInd}).subjList(subjInd) {'prob'}])
            end
        end
    end
    
    
    
    acc = res.(spaceList{spaceInd}).acc;
    nObs = res.(spaceList{spaceInd}).nObs;
%     p = res.(spaceList{spaceInd}).p;
    
    acc = sum(acc.*nObs,2)./sum(nObs,2);
    nObs = sum(nObs,2);
    acc(end+1) = sum(acc.*nObs,1)./sum(nObs,1);
    nObs(end+1) = sum(nObs,1);
    p = binocdf(acc.*nObs,nObs,0.5,'upper');
    [P,~,STATS] = signrank(acc(1:end-1),0.5,'tail','right');
    [~,pci] = binofit(round(acc.*nObs),nObs,0.1);
%     [~,pci] = binofit(acc.*nObs,nObs,0.05);
    
    accGroup(:,spaceInd) = acc(end);
    accSubj(:,spaceInd) = acc(1:end-1); clear acc
    nObsSubj(:,spaceInd) = nObs(1:end-1);
    nObsGroup(:,spaceInd) = nObs(end); clear nObs
    pSubj(:,spaceInd) = p(1:end-1);
    pGroup(:,spaceInd) = p(end); clear p
    accSubj5(:,spaceInd) = pci(1:end-1,1);
    accSubj95(:,spaceInd) = pci(1:end-1,2);
    accGroup5(1,spaceInd) = pci(end,1);
    accGroup95(1,spaceInd) = pci(end,2); clear pci
    
    signedrank(1,spaceInd) = STATS.signedrank; clear STATS
    signedrankP(1,spaceInd) = P; clear P
end
disp(spaceList)
disp(['signedrank   = ' num2str(signedrank,'%0.2f  ')])
disp(['signedrank p = ' num2str(signedrankP,'%0.4f  ')])
disp(['group accuracy  = ' num2str(accGroup*100,'%0.2f%%  ')])
disp(['binomial   p    = ' num2str(pGroup,'%0.4f   ')])

%% Group stats
switch groupStatMethod
    case 'pseudoMedian' %90%CI (90% because one-sided test)
        accPseudoMedian = pseudomedian(accSubj,1);
        accPseudoMedian5 = accPseudoMedian;
        accPseudoMedian95 = accPseudoMedian;
        delta = 0.000001;
        done = 0; i = 0;
        while ~done
            i = i+1;
            
            M = accPseudoMedian-i*delta;
            for spaceInd = 1:length(spaceList)
                [P5(spaceInd),~,~] = signrank(accSubj(:,spaceInd),M(spaceInd),'tail','right');
            end
            accPseudoMedian5(P5>=0.05) = M(P5>=0.05);
            
            M = accPseudoMedian+i*delta;
            for spaceInd = 1:length(spaceList)
                [P95(spaceInd),~,~] = signrank(accSubj(:,spaceInd),M(spaceInd),'tail','left');
            end
            accPseudoMedian95(P95>=0.05) = M(P95>=0.05);
            
            if all([P95 P5]<0.05); done=1; end
            groupE = accPseudoMedian;
            groupCI5 = accPseudoMedian5;
            groupCI95 = accPseudoMedian95;
        end
    case 'binomial'
        groupE = accGroup;
        groupCI5 = accGroup5;
        groupCI95 = accGroup95;
    otherwise
        error('X')
end

%% Pseudomedian 

f = figure('WindowStyle','docked');
y = cat(1,accSubj,groupE)';
hb = bar(y);
ylim([0 1])
hl = line(xlim,[1 1].*0.5,'color','k');
uistack(hl,'bottom')
hold on

i = length(hb);
set(hb(i),'FaceColor',[1 1 1].*0)
if plotUpperErrorBar
    heb(i) = errorbar(hb(i).XEndPoints,groupE,groupE-groupCI5,groupCI95-groupE,'.');
else
    heb(i) = errorbar(hb(i).XEndPoints,groupE,groupE-groupCI5,[],'.');
end
heb(i).Marker = 'none';
heb(i).LineWidth = lw;
heb(i).CapSize = 0;
heb(i).Color = 'r';
% heb(i).Color = [1 1 1].*0.5;


for subjInd = 1:length(hb)-1
    set(hb(subjInd),'FaceColor',[1 1 1].*0.9)
    if plotUpperErrorBar
        heb(subjInd) = errorbar(hb(subjInd).XEndPoints,accSubj(subjInd,:),accSubj(subjInd,:)-accSubj5(subjInd,:),accSubj95(subjInd,:)-accSubj(subjInd,:),'.');
    else
        heb(subjInd) = errorbar(hb(subjInd).XEndPoints,accSubj(subjInd,:),accSubj(subjInd,:)-accSubj5(subjInd,:),[],'.');
    end
    heb(subjInd).Marker = 'none';
    heb(subjInd).LineWidth = lw;
    heb(subjInd).CapSize = 0;
    heb(subjInd).Color = 'k';
end
barWidth = (hb(2).XEndPoints(1) - hb(1).XEndPoints(1)) * hb(1).BarWidth;

ax = gca;
ax.XAxis.TickLabelInterpreter = 'none';
ax.XTick = 1:length(spaceList);
ax.XTickLabel = spaceList;
ax.TickLength = [0 0];

y = 0.01;
for subjInd = 1:length(subjList)
    x = hb(subjInd).XEndPoints;
    ht(subjInd,:) = text(x-barWidth*0.05,ones(size(spaceList)).*y,subjList{subjInd},'Rotation',90,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',8);
end
subjInd = length(subjList)+1;
x = hb(subjInd).XEndPoints;
switch groupStatMethod
    case 'pseudoMedian'
        ht(subjInd,:) = text(x-barWidth*0.05,ones(size(spaceList)).*y,'median','Rotation',90,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',8,'Color','w');
    case 'binomial'
        ht(subjInd,:) = text(x-barWidth*0.05,ones(size(spaceList)).*y,'all participants','Rotation',90,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',8,'Color','w');
    otherwise
        error('x')
end


ax.YTick = 0:0.05:1;
ax.YTickLabel = cellstr(num2str(ax.YTick'*100,'%3.0f%%'));
ax.YLabel.String = {'Decoding Accuracy' '(5% lower bound of binomial distribution)'};
ax.XLabel.String = {'Decoded Response Features'};
ax.Title.String = 'Decoding Brain Response for Stimulus Orientation';

tmp = findobj(ax.Children,'type','Text');
if ischar(yLim)
    yPos_label = 0;
else
    yPos_label = yLim(1);
end
for i = 1:length(tmp)
    tmp(i).Position(2) = tmp(i).Position(2)+yPos_label-0.005;
end
ylim(yLim)
ax.YGrid = 'on';
ax.Box = 'off';

if saveFig
    filename = fullfile(pwd,mfilename);
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,'acc');
    f.Color = 'none';
    set(findobj(f.Children,'type','Axes'),'color','none')
    saveas(f,[filename '.svg']); disp([filename '.svg'])
    f.Color = 'w';
    set(findobj(f.Children,'type','Axes'),'color','w')
    saveas(f,filename); disp([filename '.fig'])
end


