function plotDecoding(res)
spaceList = fields(res)';
subjList = res.(spaceList{1}).subjList;

for spaceInd = 1:length(spaceList)
    acc = res.(spaceList{spaceInd}).acc;
    nObs = res.(spaceList{spaceInd}).nObs;
    p = res.(spaceList{spaceInd}).p;
    
    acc = sum(acc.*nObs,2)./sum(nObs,2);
    nObs = sum(nObs,2);
    p = binocdf(acc.*nObs,nObs,0.5,'upper');
    [P,~,STATS] = signrank(acc,0.5,'tail','right');
    [~,pci] = binofit(acc.*nObs,nObs,0.05);
    
    accSubj(:,spaceInd) = acc;
    nObsSubj(:,spaceInd) = nObs;
    pSubj(:,spaceInd) = p;
    accSubj5(:,spaceInd) = pci(:,1);
    accSubj95(:,spaceInd) = pci(:,2);
    
    rankGroup(1,spaceInd) = STATS.signedrank;
    pGroup(1,spaceInd) = P;
    
end
rankGroup
pGroup

%% Pseudomedian 90%CI (because one-sided test)
accPseudoMedian = pseudomedian(accSubj,1);
accPseudoMedian5 = accPseudoMedian;
accPseudoMedian95 = accPseudoMedian;
delta = 0.00001;
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
end

figure('WindowStyle','docked');
y = cat(1,accSubj,accPseudoMedian)';
hb = bar(y);
ylim([0 1])
hl = line(xlim,[1 1].*0.5,'color','k');
uistack(hl,'bottom')
hold on

i = 7;
set(hb(i),'FaceColor',[1 1 1].*0)
heb(i) = errorbar(hb(i).XEndPoints,accPseudoMedian,accPseudoMedian-accPseudoMedian5,accPseudoMedian95-accPseudoMedian,'.');
heb(i).Marker = 'none';
heb(i).CapSize = 0;
heb(i).Color = 'r';


for subjInd = 1:6
    set(hb(subjInd),'FaceColor',[1 1 1].*0.8)
    heb(subjInd) = errorbar(hb(subjInd).XEndPoints,accSubj(subjInd,:),accSubj(subjInd,:)-accSubj5(subjInd,:),accSubj95(subjInd,:)-accSubj(subjInd,:),'.');
    heb(subjInd).Marker = 'none';
    heb(subjInd).CapSize = 0;
    heb(subjInd).Color = 'k';
end
barWidth = (hb(2).XEndPoints(1) - hb(1).XEndPoints(1)) * hb(1).BarWidth;

ax = gca;
ax.XTickLabel = spaceList;
ax.TickLength = [0 0];

y = 0.01;
for subjInd = 1:length(subjList)
    x = hb(subjInd).XEndPoints;
    ht(subjInd,:) = text(x-barWidth*0.05,ones(size(spaceList)).*y,subjList{subjInd},'Rotation',90,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',8);
end
subjInd = length(subjList)+1;
x = hb(subjInd).XEndPoints;
ht(subjInd,:) = text(x-barWidth*0.05,ones(size(spaceList)).*y,'median','Rotation',90,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',8,'Color','w');

ax.YTick = 0:0.25:1;
ax.YTickLabel = cellstr(num2str(ax.YTick'*100));
ax.YLabel.String = {'Decoding Accuracy' '(%correct +/-CI)'};
ax.XLabel.String = {'Decoded Response Features'};
ax.Title.String = 'Decoding Brain Response for Stimulus Orientation';




