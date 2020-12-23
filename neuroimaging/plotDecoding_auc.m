function plotDecoding_auc(res,figOption)
groupStatMethod = 'pseudoMedian'; % 'binomial' or 'pseudoMedian'
if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 1;
    figOption.subj = 1; % 'all' or subjInd
end

plotUpperErrorBar = 0;
lw = 1;
yLim = [0 100]./100;
% yLim = 'auto';

spaceList = fields(res)';
subjList = res.(spaceList{1}).subjList;

for spaceInd = 1:length(spaceList)
    auc = res.(spaceList{spaceInd}).auc;
    nObs = res.(spaceList{spaceInd}).nObs;
    auc = auc.*nObs./mean(nObs,2); % weight within-subject accroding to number of obs in each session
    auc = mean(auc,2);
    nObs = sum(nObs,2);
    auc(end+1) = mean(auc,1);
    nObs(end+1) = sum(nObs,1);
    [P,~,STATS] = signrank(auc(1:end-1),0.5,'tail','right');

    aucGroup(:,spaceInd) = auc(end);
    aucSubj(:,spaceInd) = auc(1:end-1); clear acc
    nObsSubj(:,spaceInd) = nObs(1:end-1);
    nObsGroup(:,spaceInd) = nObs(end); clear nObs
    
    signedrank(1,spaceInd) = STATS.signedrank; clear STATS
    signedrankP(1,spaceInd) = P; clear P
end
disp(spaceList)
disp(['signedrank                  = ' num2str(signedrank,'%0.2f  ')])
disp(['signedrank p                = ' num2str(signedrankP,'%0.4f  ')])
disp(['group aera under the curve  = ' num2str(aucGroup*100,'%0.2f  ')])

%% Group stats
switch groupStatMethod
    case 'pseudoMedian' %90%CI (90% because one-sided test)
        aucPseudoMedian = pseudomedian(aucSubj,1);
        aucPseudoMedian5 = aucPseudoMedian;
        aucPseudoMedian95 = aucPseudoMedian;
        delta = 0.0001;
        done = 0; i = 0;
        while ~done
            i = i+1;

            M = aucPseudoMedian-i*delta;
            for spaceInd = 1:length(spaceList)
                [P5(spaceInd),~,~] = signrank(aucSubj(:,spaceInd),M(spaceInd),'tail','right');
            end
            aucPseudoMedian5(P5>=0.05) = M(P5>=0.05);

            M = aucPseudoMedian+i*delta;
            for spaceInd = 1:length(spaceList)
                [P95(spaceInd),~,~] = signrank(aucSubj(:,spaceInd),M(spaceInd),'tail','left');
            end
            aucPseudoMedian95(P95>=0.05) = M(P95>=0.05);

            if all([P95 P5]<0.05); done=1; end
            groupE = aucPseudoMedian;
            groupCI5 = aucPseudoMedian5;
            groupCI95 = aucPseudoMedian95;
        end
    case 'binomial'
        groupE = aucGroup;
        groupCI5 = accGroup5;
        groupCI95 = accGroup95;
    otherwise
        error('X')
end

%% Pseudomedian

f = figure('WindowStyle','docked');
y = cat(1,aucSubj,groupE)';
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
%     if plotUpperErrorBar
%         heb(subjInd) = errorbar(hb(subjInd).XEndPoints,aucSubj(subjInd,:),aucSubj(subjInd,:)-accSubj5(subjInd,:),accSubj95(subjInd,:)-aucSubj(subjInd,:),'.');
%     else
%         heb(subjInd) = errorbar(hb(subjInd).XEndPoints,aucSubj(subjInd,:),aucSubj(subjInd,:)-accSubj5(subjInd,:),[],'.');
%     end
%     heb(subjInd).Marker = 'none';
%     heb(subjInd).LineWidth = lw;
%     heb(subjInd).CapSize = 0;
%     heb(subjInd).Color = 'k';
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


% ax.YTick = 0:0.05:1;
% ax.YTickLabel = cellstr(num2str(ax.YTick','%0.1f'));
ax.YLabel.String = {'Area Under the ROC'};
ax.XLabel.String = {'Decoded Response Features'};
ax.Title.String = 'Decoding Brain Response for Stimulus Orientation';

% tmp = findobj(ax.Children,'type','Text');
% if ischar(yLim)
%     yPos_label = 0;
% else
%     yPos_label = yLim(1);
% end
% for i = 1:length(tmp)
%     tmp(i).Position(2) = tmp(i).Position(2)+yPos_label-0.005;
% end
ylim(yLim)
ax.YGrid = 'on';
ax.Box = 'off';

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
    writeFig(f,mfilename,'auc')
end
