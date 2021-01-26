function plotDecoding_acc(res,figOption,verbose,accPerm)
groupStatMethod = 'binomial'; % 'binomial' or 'pseudoMedian'
if ~exist('verbose','var')
    verbose = 1;
end
if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end

plotUpperErrorBar = 1;
lw = 1;
yLim = [0 100]./100;
% yLim = 'auto';

spaceList = fields(res)';
subjList = res.(spaceList{1}).subjList;

for spaceInd = 1:length(spaceList)
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
t = table(...
    num2str(accGroup'*100,'%0.2f%%'),...
    num2str(pGroup','%0.4f'),...
    num2str(signedrank','%0.2f'),...
    num2str(signedrankP','%0.4f'));
t.Properties.VariableNames = {'acc' 'bino p' 'signed rank' 'signed rank one-tailed p'};
t.Properties.RowNames = spaceList;
disp(t)

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
% ylim([0 1])
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
barWidthAx = (hb(2).XEndPoints(1) - hb(1).XEndPoints(1)) * hb(1).BarWidth;
ylim(yLim)
ax = gca;
ax.XAxis.TickLabelInterpreter = 'none';
ax.XTick = 1:length(spaceList);
ax.XTickLabel = spaceList;
ax.TickLength = [0 0];

y = 0.01;
for subjInd = 1:length(subjList)
    x = hb(subjInd).XEndPoints;
    ht(subjInd,:) = text(x-barWidthAx*0.05,ones(size(spaceList)).*y,subjList{subjInd},'Rotation',90,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',8);
end
subjInd = length(subjList)+1;
x = hb(subjInd).XEndPoints;
switch groupStatMethod
    case 'pseudoMedian'
        ht(subjInd,:) = text(x-barWidthAx*0.05,ones(size(spaceList)).*y,'median','Rotation',90,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',8,'Color','w');
    case 'binomial'
        ht(subjInd,:) = text(x-barWidthAx*0.05,ones(size(spaceList)).*y,'all participants','Rotation',90,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',8,'Color','w');
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

ax.YGrid = 'on';
ax.Box = 'off';
drawnow

if exist('accPerm','var')
    pos = cell2mat(get(hb,'XEndPoints'));
    deltaGroup = diff(sort([pos(1,:) pos(end,:)])); deltaGroup = mode(deltaGroup(2:2:end));
    ax.XLim = [ax.XLim(1) pos(end)+deltaGroup];
    drawnow
    
    barWidth = mode(diff(pos,[],1)) .* mode(cell2mat(get(hb,'BarWidth')));
    barWidthPos = barWidth/diff(ax.XLim)*ax.Position(3);
    deltaGroupPos = deltaGroup/diff(ax.XLim)*ax.Position(3);
    
    pPerm = nan(size(spaceList));
    for spaceInd = 1:length(spaceList)
        nPerm = size(accPerm.(spaceList{spaceInd}),1);
        nObs = sum(res.(spaceList{spaceInd}).nObs(:));
        threshPerm = quantile(accPerm.(spaceList{spaceInd}),[0.05 1-0.05]).*nObs;
        threshBino = binoinv([0.05 1-0.05],nObs,0.5);
%         plot(x,thresh,'r_')
        pPerm(spaceInd) = sum(accPerm.(spaceList{spaceInd})>accGroup(spaceInd))/nPerm;
        
        hit = round(accPerm.(spaceList{spaceInd}).*nObs);
        xRatio = ( hb(7).XEndPoints(spaceInd) + barWidth(spaceInd)/2 - ax.XLim(1) ) / diff(ax.XLim);
        xRatio = xRatio*ax.Position(3);
        axHist(spaceInd) = axes; drawnow
        axHist(spaceInd).XTick = []; drawnow
        height = deltaGroupPos-barWidthPos(spaceInd);
        axHist(spaceInd).Position([1 3]) = [ax.Position(1)+xRatio height] + [0 -0.25].*height;
        axHist(spaceInd).Position([2 4]) = ax.Position([2 4]);
        plotHitsAndBino(hit,nObs); drawnow
        ylim([0 max([axHist(spaceInd).Children(1).YData(:); axHist(spaceInd).Children(2).YData(:)])]); drawnow
        axHist(spaceInd).XDir = 'reverse';
        view(axHist(spaceInd),[90 90]);
        axHist(spaceInd).XLim = ax.YLim.*nObs; drawnow
        axHist(spaceInd).Box = 'off';
        axHist(spaceInd).Color = 'none';
        hLine = findobj(axHist(spaceInd).Children,'type','line');
        hBar = findobj(axHist(spaceInd).Children,'type','Patch');
        ind = hLine.XData > threshBino(2);
        hLine.YData(ind) = [];
        hLine.XData(ind) = [];
        ind = hLine.XData < threshBino(1);
        hLine.YData(ind) = [];
        hLine.XData(ind) = [];
        ind = mean(hBar.XData,1) > threshPerm(2);
        hBar.YData(:,ind) = [];
        hBar.XData(:,ind) = [];
        ind = mean(hBar.XData,1) < threshPerm(1);
        hBar.YData(:,ind) = [];
        hBar.XData(:,ind) = [];
        hLine.YData(end+1) = 0;
        hLine.XData(end+1) = hLine.XData(end);
        hLine.YData(end+1) = 0;
        hLine.XData(end+1) = hLine.XData(1);
        hLine.YData(end+1) = hLine.YData(1);
        hLine.XData(end+1) = hLine.XData(1);
        hBar.XData(1:2,1) = mean(hBar.XData(:,1));
        hBar.XData(3:4,end) = mean(hBar.XData(:,end));
        axHist(spaceInd).XAxis.Visible = 'off';
        axHist(spaceInd).YAxis.Visible = 'off';
%         plot([0 threshBino],[0 0],'-r')
%         plot(axHist(spaceInd).XLim,[0 0],'w')
        drawnow
    end
    for spaceInd = 1:length(spaceList)
        xRatio = ( hb(7).XEndPoints(spaceInd) + barWidth(spaceInd)/2 - ax.XLim(1) ) / diff(ax.XLim);
        xRatio = xRatio*ax.Position(3);
        height = deltaGroupPos-barWidthPos(spaceInd);
        axHist(spaceInd).Position(1) = ax.Position(1)+xRatio + 0.03*height;
    end
end

switch groupStatMethod
    case 'pseudoMedian'
        tmpStr1 = 'Pseudo-median';
        tmpStr2 = 'CI (5% lower bound)';
    case 'binomial'
        tmpStr1 = 'Mean';
        tmpStr2 = 'Binomial CI (5% lower bound)';
    otherwise
        error('X')
end
legend(ax.Children([end-1 end-7 end-8]),char({'Participants' tmpStr1 tmpStr2}),'Box','off');


if figOption.save
    writeFig(f,mfilename,'acc',verbose)
end
