function f = plotAllDecoding(p,res,info)

nBoot = 2^13;
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
yInfo = 'respFeat x condPair x subj';
% y = y([2 3 1],:,:);
% respFeatList = respFeatList([2 3 1]);

y = cat(2,y(:,1,:),mean(y(:,2:3,:),2));
% y = cat(1,y(1,:,:),mean(y(2:3,:,:),1));
condPairList{2} = 'gratVSplaid';
condPairList(3) = [];

yAv = mean(y,3);
[yBoot_90CI,~] = bootThat(y,nBoot);
yErLow = yAv-yBoot_90CI(:,:,1);
yErHigh = yBoot_90CI(:,:,2)-yAv;
% yEr = std(y,[],3);

f = figure('WindowStyle','docked');
hBar = bar(yAv); hold on
for i = 1:length(condPairList)
% for i = 1:length(respFeatList)
    x = hBar(i).XEndPoints;
    hErBar = errorbar(x,yAv(:,i),yErLow(:,i),yErHigh(:,i));
    hErBar.LineStyle = 'none';
    hErBar.Color = 'k';
    hErBar.CapSize = 0;
end
ylim([0 1]);

ax = gca;
% ax.XTickLabel = condPairList;
ax.XTickLabel = respFeatList;
hP = plot(xlim,[0.5 0.5],'k');
uistack(hP,'bottom')
% legend(hBar,respFeatList)
legend(hBar,condPairList)

%% Save
if p.figOption.save
    fullfilename = fullfile(p.figOption.finalDir,'decodingSummary');
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
sz = size(y);
sz(1) = nBoot;
yBoot = nan(sz);
for bootInd = 1:nBoot
    i = nan(nSubj,1);
    for drawInd = 1:nSubj
        i(drawInd) = randperm(nSubj,1);
    end
    yBoot(bootInd,:) = mean(y(i,:),1);
end

sz = size(yBoot);
yBoot = permute(yBoot,[2:length(sz) 1]);
sz = size(yBoot);
yBoot_90CI = prctile(yBoot,[5 95],length(sz));
yBoot_95CI = prctile(yBoot,[2.5 97.5],length(sz));


