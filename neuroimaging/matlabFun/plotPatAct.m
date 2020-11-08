function [h,hPol] = plotPatAct(svm,h,hPol,plotIt)
if ~exist('plotIt','var')
    plotIt = 1;
end
if ~exist('h','var') || isemptry(h)
    h = nan(1,6);
end
if ~exist('hPol','var') || isemptry(hPol)
    hPol = nan(1,length(svm.p.trainInfo));
end
titleStr = svm.p.dataFileOut;
tmpLength = length(titleStr);
titleStr = {titleStr(1:round(tmpLength/2)); titleStr(round(tmpLength/2)+1:end)};




%% Plot acc
% Initiate plot
if isnan(h(1))
    if plotIt
        h(1) = figure('windowStyle','docked');
    else
        h(1) = figure('visible','off');
    end
    colormap jet
else
    if plotIt
        figure(h(1)); hold off;
    else
        figure(h(1),'visible','off'); hold off;
    end
end

% Compute accuracies
acc = squeeze(nanmean(nanmean(svm.r.hitRate,1),4));

% Plot it
imagesc(acc,[0 1]); colorbar
xlabel('test')
ylabel('train')
Tick = 1:length(acc);
set(gca,'YTick',Tick)
set(gca,'YTickLabel',svm.r.crossInfo_dim1Train)
set(gca,'XTick',Tick)
set(gca,'XTickLabel',svm.r.crossInfo_dim2Test)
set(gca, 'XAxisLocation', 'top')
set(gca, 'TickDir', 'out')
title([titleStr; {'Accuracies'}],'interpret','none')



%% Plot acc thresholded
% Initiate plot
if isnan(h(2))
    if plotIt
        h(2) = figure('windowStyle','docked');
    else
        h(2) = figure('visible','off');
    end
    colormap jet
else
    if plotIt
        figure(h(2)); hold off;
    else
        figure(h(2),'visible','off'); hold off;
    end
end

% Threshold accuracies
negThresh = binoinv(0.05,svm.p.nObs,0.5)/svm.p.nObs;
posThresh = binoinv(1-0.05,svm.p.nObs,0.5)/svm.p.nObs;
acc(acc>negThresh & acc<posThresh) = nan;

% Plot it
imagesc(acc,[0 1]); colorbar
xlabel('test')
ylabel('train')
Tick = 1:length(acc);
set(gca,'YTick',Tick)
set(gca,'YTickLabel',svm.r.crossInfo_dim1Train)
set(gca,'XTick',Tick)
set(gca,'XTickLabel',svm.r.crossInfo_dim2Test)
set(gca, 'XAxisLocation', 'top')
set(gca, 'TickDir', 'out')
title([titleStr; {'Thresholded Accuracies'}],'interpret','none')





%% Plot decision value
% Compute d
d = squeeze(mean(svm.r.crossInfo.d,1));
d = cat(4,d(:,:,1:end/2),d(:,:,end/2+1:end));
%Compute magnitude of the difference
dDiff = -mean(diff(d,[],4),3);
%Compute statistical difference
for i = 1:size(d,1)
    for ii = 1:size(d,2)
        dTmp = squeeze(d(i,ii,:,:));
        [H,P,CI,STATS] = ttest2(dTmp(:,1),dTmp(:,2));
        dT(i,ii) = STATS.tstat;
        dP(i,ii) = P;
    end
end
dDiff_thresh = dDiff;
dDiff_thresh(dP>0.05) = nan;

% Initiate plot
if isnan(h(4))
    if plotIt
        h(4) = figure('windowStyle','docked');
    else
        h(4) = figure('visible','off');
    end
    colormap jet
else
    if plotIt
        figure(h(4)); hold off;
    else
        figure(h(4),'visible','off'); hold off;
    end
end
% Plot dDiff
imagesc(dDiff_thresh); colorbar
xlabel('test')
ylabel('train')
Tick = 1:length(acc);
set(gca,'YTick',Tick)
set(gca,'YTickLabel',svm.r.crossInfo_dim1Train)
set(gca,'XTick',Tick)
set(gca,'XTickLabel',svm.r.crossInfo_dim2Test)
set(gca, 'XAxisLocation', 'top')
set(gca, 'TickDir', 'out')
title([titleStr; {'Diff in decision; thresholded'}],'interpret','none')


% Initiate plot
if isnan(h(5))
    if plotIt
        h(5) = figure('windowStyle','docked');
    else
        h(5) = figure('visible','off');
    end
    colormap jet
else
    if plotIt
        figure(h(5)); hold off;
    else
        figure(h(5),'visible','off'); hold off;
    end
end
% Plot dDiff thresholded
imagesc(dDiff); colorbar
xlabel('test')
ylabel('train')
Tick = 1:length(acc);
set(gca,'YTick',Tick)
set(gca,'YTickLabel',svm.r.crossInfo_dim1Train)
set(gca,'XTick',Tick)
set(gca,'XTickLabel',svm.r.crossInfo_dim2Test)
set(gca, 'XAxisLocation', 'top')
set(gca, 'TickDir', 'out')
title([titleStr; {'Diff in decision'}],'interpret','none')




% Initiate plot
if isnan(h(6))
    if plotIt
        h(6) = figure('windowStyle','docked');
    else
        h(6) = figure('visible','off');
    end
    colormap jet
else
    if plotIt
        figure(h(6)); hold off;
    else
        figure(h(6),'visible','off'); hold off;
    end
end
% Plot dT
imagesc(dT); colorbar
xlabel('test')
ylabel('train')
Tick = 1:length(acc);
set(gca,'YTick',Tick)
set(gca,'YTickLabel',svm.r.crossInfo_dim1Train)
set(gca,'XTick',Tick)
set(gca,'XTickLabel',svm.r.crossInfo_dim2Test)
set(gca, 'XAxisLocation', 'top')
set(gca, 'TickDir', 'out')
title([titleStr; {'T stats of difference in decision'}],'interpret','none')






% %% Plot filtered responses
% for ind = 1:length(svm.p.trainInfo)
%     curPatAct = squeeze(svm.r.crossInfo.f(:,ind,1,:));
%     curPatAct = mean(curPatAct,1);
%     curPatAct = reshape(curPatAct,[length(curPatAct)/2 2]);
%     patAct(:,:,ind) = curPatAct;
%     
%     % Polar scatter
%     if isnan(hPol(ind))
%         if plotIt
%             hPol(ind) = figure('windowStyle','docked');
%         else
%             hPol(ind) = figure('visible','off');
%         end
%         colormap jet
%     else
%         if plotIt
%             figure(hPol(ind)); hold off;
%         else
%             figure(hPol(ind),'visible','off'); hold off;
%         end
%     end
%     
%     scatter(real(curPatAct(:,1)),imag(curPatAct(:,1)),'k'); hold on
%     scatter(real(curPatAct(:,2)),imag(curPatAct(:,2)),'r')
%     plot(get(gca,'xlim'),[0 0],'k')
%     plot([0 0],get(gca,'ylim'),'k')
%     
%     
%     xlabel('real')
%     ylabel('imag')
%     legend({'ori1','ori2'})
%     title([titleStr; {['Trained on ' svm.p.trainInfo{ind}]}],'interpret','none')
% end
% 
% % Bar plot of amp
% if isnan(h(6))
%     if plotIt
%         h(6) = figure('windowStyle','docked');
%     else
%         h(6) = figure('visible','off');
%     end
% else
%     if plotIt
%         figure(h(6)); hold off;
%     else
%         figure(h(6),'visible','off'); hold off;
%     end
% end
% 
% curPatAct = abs(patAct);
% curMean = squeeze(mean(curPatAct,1))';
% curStd = squeeze(std(curPatAct,[],1) / sqrt(size(curPatAct,1)))';
% bar(curMean); hold on
% curX = [get(gca,'XTick')'-0.15 get(gca,'XTick')'+0.15];
% errorbar(curX,curMean,curStd,'.')
% set(gca,'XTick',1:length(svm.p.trainInfo));set(gca,'XTickLabel',svm.p.trainInfo)
% ylabel('amp (+-sem)')
% title([titleStr; {['']}],'interpret','none')
% 
% 
% % Bar plot of delay
% if isnan(h(7))
%     if plotIt
%         h(7) = figure('windowStyle','docked');
%     else
%         h(7) = figure('visible','off');
%     end
%     colormap jet
% else
%     if plotIt
%         figure(h(7)); hold off;
%     else
%         figure(h(7),'visible','off'); hold off;
%     end
% end
% 
% curPatAct = angle(patAct);
% curMean = squeeze(circ_mean(curPatAct,[],1))'/(pi/2)*6;
% curStd = squeeze(circ_std(curPatAct,[],1))'/(pi/2)*6 / sqrt(size(curPatAct,1));
% bar(curMean); hold on
% curX = [get(gca,'XTick')'-0.15 get(gca,'XTick')'+0.15];
% errorbar(curX,curMean,curStd,'.')
% set(gca,'XTick',1:length(svm.p.trainInfo));set(gca,'XTickLabel',svm.p.trainInfo)
% ylabel('delay (sec +- sem)')
% title([titleStr; {['']}],'interpret','none')









