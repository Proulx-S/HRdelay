function f = plotAllDecoding2(p,res,info)
p.figOption.save = 0;
verbose  = 0;

nBoot = p.boot.n;
metric = 'auc';
condPairList = info.condPairList;
respFeatList = info.respFeatList;
n = size(res{1,1}.subj.(metric),1);
y = nan(length(condPairList),length(respFeatList),n);
for respFeatInd = 1:length(respFeatList)
    for condPairInd = 1:length(condPairList)
        y(condPairInd,respFeatInd,:) = res{condPairInd,respFeatInd}.subj.(metric);
    end
end

% yInfo = 'condPair x respFeat x subj';
y = permute(y,[2 1 3]);
condPairList = condPairList';
respFeatList = respFeatList';
yInfo = 'respFeat x condPair x subj';
% y = y([2 3 1],:,:);
% respFeatList = respFeatList([2 3 1]);

y2 = mean(y,2);
% y = cat(1,y(1,:,:),mean(y(2:3,:,:),1));
condPairList2 = condPairList;
% condPairList2{2} = 'gratVSplaid';
% condPairList2(3) = [];

yAv = mean(y2,3);
[yBoot_90CI,~] = bootThat(y2,nBoot);
yErLow = yAv-yBoot_90CI(:,:,1);
yErHigh = yBoot_90CI(:,:,2)-yAv;
% yEr = std(y,[],3);

f = figure('WindowStyle','docked');
hBar = bar(yAv); hold on
for i = 1%:length(condPairList2)
% for i = 1:length(respFeatList)
    x = hBar.XEndPoints;
    hErBar = errorbar(x,yAv(:,i),yErLow(:,i),yErHigh(:,i));
    hErBar.LineStyle = 'none';
    hErBar.Color = 'k';
    hErBar.CapSize = 0;
end
ylim([0 1]);
ax = gca;






if isfield(res{1}.subj,'aucP')
    tmp = reshape([res{:}],size(res));
    tmp = reshape([tmp.subj],size(tmp));
    tmp = reshape(permute([tmp.aucP],[2 3 1]),[size(tmp) size(tmp(1).aucP,1)]);
    aucP = permute(tmp,[3 1 2]);
    aucP = mean(aucP,2);
    info.condPairList = {'all'};
    
    fTmp = figure('windowstyle','docked','visible',verbose);
    binWidth = 5;
    t = tiledlayout(3,1,'TileIndexing','columnmajor');
    binN = cell(size(aucP,[2 3]));
    binEdges = cell(size(aucP,[2 3]));
    for respFeatInd = 1:size(aucP,3)
        axTmp(respFeatInd) = nexttile;
        for condPairInd = 1:size(aucP,2)
            hHist(condPairInd) = histogram(aucP(:,condPairInd,respFeatInd),linspace(prctile(aucP(:,condPairInd,respFeatInd),5),prctile(aucP(:,condPairInd,respFeatInd),95),100/binWidth-2*5/binWidth),'Normalization','pdf'); hold on
            [binN{condPairInd,respFeatInd},binEdges{condPairInd,respFeatInd}] = histcounts(aucP(:,condPairInd,respFeatInd),linspace(prctile(aucP(:,condPairInd,respFeatInd),5),prctile(aucP(:,condPairInd,respFeatInd),95),100/binWidth-2*5/binWidth),'Normalization','pdf');
            hHist(condPairInd).EdgeColor = 'none';
        end
        xlim([0 1])
        grid on
        hLeg(respFeatInd) = legend(info.condPairList);
        title(respFeatList{respFeatInd})
    end
    linkaxes(axTmp,'x');
    xticklabels(axTmp(1:2),{})
    xlabel(t,'auc')
    title(t,'Permutation Test')
    t.TileSpacing = 'compact';
    set(hLeg(2:3),'Visible','off')
    hLeg(1).Box = 'off';
    hLeg(1).Location = 'northwest';
    hLeg(1).Position(1) = axTmp(1).Position(1);
    hLeg(1).Position(3) = 0.2;
    % sum(hHist(1).Values).*hHist(1).BinWidth
    
    
    
    
    
    figure(f)
    BarWidth = 0.5;
    scaleDown = -0.025;
    set(hBar,'BarWidth',BarWidth)
    hBarP = cell(length(info.condPairList),length(info.respFeatList));
    for respFeatInd = 1:length(info.respFeatList)
        barCent = hBar.XData(respFeatInd);
        base = barCent - hBar.BarWidth/2;
        width = mode(diff(binEdges{respFeatInd}(2:end)));
        n = binN{respFeatInd}*scaleDown;
        cent = binEdges{respFeatInd}(2:end)-width/2;
        hBarP{respFeatInd} = barh(cent,n,'hist');
        hBarP{respFeatInd}.EdgeColor = 'none';
        hBarP{respFeatInd}.FaceColor = hBar.FaceColor;
        hBarP{respFeatInd}.Vertices(:,1) = hBarP{respFeatInd}.Vertices(:,1) + base;
        uistack(hBarP{respFeatInd},'bottom')
    end
    
    
    
    axLim = axis;
    hTex1 = text(axLim(1),axLim(4),[num2str(size(res{1}.subj.aucP,1)) ' permutations/bootstraps'],'VerticalAlignment','top');
    hTex2 = text(axLim(2),axLim(4),['n=' num2str(length(res{1}.subj.subjList))],'VerticalAlignment','top','horizontalAlignment','right');
end


% x = nan(size(y));
% for respFeatInd = 1:length(respFeatList)
%     respFeatList{respFeatInd}
%     for condInd = 1:length(condPairList)
%         condPairList{condInd}
%         switch condInd
%             case 1
%                 barInd = 1;
%                 offset = 0;
%             case 2
%                 barInd = 2;
%                 offset = -abs(delta)/2;
%             case 3
%                 barInd = 2;
%                 offset = abs(delta)/2;
%             otherwise
%                 error('X')
%         end
%         x(respFeatInd,condInd,:) = hBar(barInd).XEndPoints(respFeatInd) + offset;
%     end
% end
xTmp = 1:3;%permute(x,[3 1 2]);
yTmp = permute(mean(y,2),[3 1 2]);
Mrkr = 'osd^v><ph';
for subjInd = 1:size(y,3)
    hScat(subjInd) = scatter(xTmp,yTmp(subjInd,:),Mrkr(subjInd));
end
% alpha(hScat,0.3)
set(hScat,'MarkerEdgeColor',[1 1 1].*0.4)
set(hBar,'FaceColor',[1 1 1].*0.6)
set(hBar,'EdgeColor','k')
set([hBarP{:}],'FaceColor',hBar.FaceColor)

ax = f.Children;
% ax.XTickLabel = condPairList;
ax.XTickLabel = respFeatList;
hP = plot(xlim,[0.5 0.5],'k');
uistack(hP,'bottom')
% legend(hBar,respFeatList)
% legend(hBar,condPairList2,'Box','off')


%% Save
if p.figOption.save
    fullfilename = fullfile(p.figOption.finalDir,'decodingSummary2');
    curF = f;
    curF.Color = 'none';
    set(findobj(curF.Children,'type','Axes'),'color','none')
    curFile = fullfilename;
    curExt = 'svg';
    saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
    curF.Color = 'w';
    curExt = 'fig';
    saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
    curExt = 'jpg';
    saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
end

function [yBoot_90CI,yBoot_95CI] = bootThat(y,nBoot)
nSubj = size(y,3);
sz = size(y);
y = permute(y,[length(sz) 1:(length(sz)-1)]);
sz = sz([length(sz) 1:(length(sz)-1)]);
sz(1) = nBoot;
yBoot = nan(sz);
for bootInd = 1:nBoot
    i = nan(nSubj,1);
    for drawInd = 1:nSubj
        i(drawInd) = randperm(nSubj,1);
    end
    yBoot(bootInd,:) = mean(y(i,:),1);
end

yBoot = permute(yBoot,[2:length(sz) 1]);
sz = sz([2:length(sz) 1]);
yBoot_90CI = prctile(yBoot,[5 95],length(sz));
yBoot_95CI = prctile(yBoot,[2.5 97.5],length(sz));


