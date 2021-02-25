function groupAna_sinResp(p,figOption,verbose)
if ~exist('verbose','var')
    verbose = 1;
end
if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end
colors = [  0         0.4470    0.7410
            0.8500    0.3250    0.0980
            0.9290    0.6940    0.1250];
lw = 1.5;

%% Define paths
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'd';
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp


%% Load data
dAll = cell(size(subjList,2),1);
for subjInd = 1:size(subjList,2)
    subj = subjList{subjInd};
    curFile = fullfile(funPath,inDir,[subj '.mat']);
    if verbose; disp(['loading: ' curFile]); end
    load(curFile,'res')
    dAll{subjInd} = res; clear res
end
d = dAll; clear dAll

x.sin = nan(length(subjList),3,1,2);
x.hr = nan(length(subjList),3,12,2);
for subjInd = 1:length(subjList)
    for sessInd = 1:2
        %% Average voxels
        % Between-session feature selection
        sess = ['sess' num2str(sessInd)];
        sessFeat = ['sess' num2str(~(sessInd-1)+1)];
        ind = true(size(d{subjInd}.(sess).sin,1),1);
        % activated voxels
        ind = ind & d{subjInd}.(sessFeat).featSel.F.act.p < p.act.threshVal;
        % non vein activated voxels
        veinMap = mean(d{subjInd}.(sessFeat).featSel.vein.map(:,:),2);
        ind = ind & veinMap<prctile(veinMap(ind),100-p.vein.percentile);
        
        xSin = mean(d{subjInd}.(sess).sin(ind,:,:,:,:,:),1);
        xHr = mean(d{subjInd}.(sess).hr(ind,:,:,:,:,:),1);
        
        %% Average runs
        x.sin(subjInd,:,:,sessInd) = mean(xSin,4);
        x.hr(subjInd,:,:,sessInd) = mean(xHr,4);
    end
end
x.info = 'subj x cond[grat1,grat2,plaid] x time x sess';

error('not finished')
%% Plot results
% Cartesian space
xDataNorm = xData - mean(xData(:,1:3),2) + mean(mean(xData(:,1:3),2),1);
xDataNormMean = mean(xDataNorm,1);
xDataNorm_theta = angle(xDataNorm);
xDataNorm_rho = abs(xDataNorm);
xDataNormMean_theta = angle(xDataNormMean);
xDataNormMean_rho = abs(xDataNormMean);
fGroup(1) = figure('WindowStyle','docked');
for condInd = 1:3
    h(condInd) = polarplot(xDataNorm_theta(:,condInd),xDataNorm_rho(:,condInd),'o'); hold on
    h(condInd).Color = colors(condInd,:);
    hM(condInd) = polarplot(xDataNormMean_theta(:,condInd),xDataNormMean_rho(:,condInd),'o'); hold on
    hM(condInd).MarkerFaceColor = colors(condInd,:);
    hM(condInd).MarkerEdgeColor = 'w';
end
uistack(hM,'top')
hl = legend(char({'ori1' 'ori2' 'plaid'}),'Location','east');
hl.Color = 'none';
hl.Box = 'off';
title('Group Response (subject effect removed)')
ax = gca;
ax.ThetaTickLabel = 12-ax.ThetaTick(1:end)/360*12;
ax.ThetaTickLabel(1,:) = '0 ';
ax.ThetaAxis.Label.String = {'delay' '(sec)'};
ax.ThetaAxis.Label.Rotation = 0;
ax.ThetaAxis.Label.HorizontalAlignment = 'left';
ax.RAxis.Label.String = 'amp (%BOLD)';
ax.RAxis.Label.Rotation = 80;

if figOption.save
    filename = fullfile(pwd,mfilename);
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,'cartesian');
    fGroup(1).Color = 'none';
    set(findobj(fGroup(1).Children,'type','Axes'),'color','none')
    set(findobj(fGroup(1).Children,'type','PolarAxes'),'color','none')
    curFile = filename;
    curExt = 'svg';
    saveas(fGroup(1),[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
    fGroup(1).Color = 'w';
    curExt = 'fig';
    saveas(fGroup(1),[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
    curExt = 'jpg';
    saveas(fGroup(1),[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
end

% Polar space
thetaRef = angle(mean(xData(:)));
xData_theta = wrapToPi(angle(xData) - thetaRef);
xData_theta = -xData_theta./pi*6;
thetaRef = -thetaRef./pi*6;
xData_theta = xData_theta + thetaRef; clear thetaRef
xData_thetaMean = mean(xData_theta,1);
xData_rho = abs(xData);
xData_rhoMean = mean(xData_rho,1);

fGroup(2) = figure('WindowStyle','docked');
[ax, pos] = tight_subplot(1,2);
spInd = 1;
axes(ax(spInd));
hb = bar(xData_rhoMean(:,[4 3])); hold on
hb1 = bar(hb.XData(1),hb.YData(1),'FaceColor','k');
hb2 = bar(hb.XData(2),hb.YData(2),'FaceColor',colors(3,:));
hp = plot(hb.XData,xData_rho(:,[4 3])','color',[1 1 1].*0.5);
set(hp,'LineWidth',lw);
delete(hb)
xlim([0.25 2.75])
ylabel({'Response amplitude' '(%BOLD)'})
ax(spInd).XTickLabel = {'ori' 'plaid'};
ylim([0 max(xData_rho(:)).*1.1])
ax(spInd).PlotBoxAspectRatio = [0.2 1 1];
box off
spInd = 2;
axes(ax(spInd));
hb = bar(xData_rhoMean(:,[1 2])); hold on
hb1 = bar(hb.XData(1),hb.YData(1),'FaceColor',colors(1,:));
hb2 = bar(hb.XData(2),hb.YData(2),'FaceColor',colors(2,:));
hp = plot(hb.XData,xData_rho(:,[1 2])','color',[1 1 1].*0.5);
set(hp,'LineWidth',lw);
delete(hb)
xlim([0.25 2.75])
ylabel({'Response amplitude' '(%BOLD)'})
ax(spInd).XTickLabel = {'ori1' 'ori2'};
ylim([0 max(xData_rho(:)).*1.1])
ax(spInd).PlotBoxAspectRatio = [0.2 1 1];
box off

if figOption.save
    filename = fullfile(pwd,mfilename);
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,'amp');
    fGroup(2).Color = 'none';
    set(findobj(fGroup(2).Children,'type','Axes'),'color','none')
    curFile = filename;
    curExt = 'svg';
    saveas(fGroup(2),[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
    fGroup(2).Color = 'w';
    curExt = 'fig';
    saveas(fGroup(2),[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
    curExt = 'jpg';
    saveas(fGroup(2),[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
end

fGroup(3) = figure('WindowStyle','docked');
[ax, pos] = tight_subplot(2, 1,[0],[0.1 0],[0.1 0]);
spInd = 1;
axes(ax(spInd));
hb = barh(xData_thetaMean(:,[4 3])); hold on
hb1 = barh(hb.XData(1),hb.YData(1),'FaceColor','k');
hb2 = barh(hb.XData(2),hb.YData(2),'FaceColor',colors(3,:));
hp = plot(xData_theta(:,[4 3])',hb.XData,'color',[1 1 1].*0.5);
set(hp,'LineWidth',lw);
delete(hb)
ylim([0.25 2.75]);
xlabel({'Response Delay' '(sec)'})
ax(spInd).YTickLabel = {'ori' 'plaid'};
xlim([0 max(xData_theta(:)).*1.1])
ax(spInd).PlotBoxAspectRatio = [1 0.2 1];
box off
hTxt = text(0,hb2.XData+hb2.BarWidth/2,['diff=' num2str(diff(xData_thetaMean(:,[4 3])),'%0.3fsec')],'VerticalAlignment','bottom');

spInd = 2;
axes(ax(spInd));
hb = barh(xData_thetaMean(:,[1 2])); hold on
hb1 = barh(hb.XData(1),hb.YData(1),'FaceColor',colors(1,:));
hb2 = barh(hb.XData(2),hb.YData(2),'FaceColor',colors(2,:));
hp = plot(xData_theta(:,[1 2])',hb.XData,'color',[1 1 1].*0.5);
set(hp,'LineWidth',lw);
delete(hb)
ylim([0.25 2.75]);
xlabel({'Response Delay' '(sec)'})
ax(spInd).YTickLabel = {'ori1' 'ori2'};
xlim([0 max(xData_theta(:)).*1.1])
ax(spInd).PlotBoxAspectRatio = [1 0.2 1];
box off
hTxt = text(0,hb2.XData+hb2.BarWidth/2,['diff=' num2str(diff(xData_thetaMean(:,[1 2])),'%0.3fsec')],'VerticalAlignment','bottom');


disp('***')
disp(['delay (plaid-ori) = ' num2str(diff(xData_thetaMean(:,[4 3])),'%0.3fs')])
disp(['delay (ori2-ori1) = ' num2str(diff(xData_thetaMean(:,[1 2])),'%0.3fs')])
disp(['amp (plaid-ori) = ' num2str(diff(xData_rhoMean(:,[4 3])),'%0.3f%%BOLD')])
disp(['amp (ori2-ori1) = ' num2str(diff(xData_rhoMean(:,[1 2])),'%0.3f%%BOLD')])
disp('***')

if figOption.save
    filename = fullfile(pwd,mfilename);
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,'delay');
    fGroup(3).Color = 'none';
    set(findobj(fGroup(3).Children,'type','Axes'),'color','none')
    curFile = filename;
    curExt = 'svg';
    saveas(fGroup(3),[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
    fGroup(3).Color = 'w';
    curExt = 'fig';
    saveas(fGroup(3),[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
    curExt = 'jpg';
    saveas(fGroup(3),[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
end

%% Stats
xDataMean = mean(xData(:));

amp = abs(xData);
delay = angle(xData);
delay = wrapToPi(delay - angle(xDataMean));
% delay = -(delay + angle(Xmean))/pi*6;


% On cartesian space (conservative)
disp('---------------')
disp('Cartesian Space')
disp('---------------')
disp('Ori vs Plaid:')
X = cat(1,xData(:,4),xData(:,3));
stats = T2Hot2d([real(X) imag(X)],0.05);
disp('Hotelling''s T^2 multivariate test')
disp([' T^2=' num2str(stats.T2,'%0.2f')]);
disp([' p=' num2str(stats.P,'%0.2f')]);
disp('Ori1 vs Ori2:')
X = cat(1,xData(:,1),xData(:,2));
stats = T2Hot2d([real(X) imag(X)],0.05);
disp(' Hotelling''s T^2 multivariate test')
disp([' T^2=' num2str(stats.T2,'%0.2f') '; p=' num2str(stats.P,'%0.2f')]);

% Amplitude in polar space
disp('---------------')
disp('Polar Amplitude')
disp('---------------')
disp('Ori vs Plaid (one-tailed):')
x = amp(:,4); y = amp(:,3);
[H,P,CI,STATS] = ttest(x,y,'tail','right');
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,H,STATS] = signrank(x,y,'tail','right');
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);
disp('Ori1 vs Ori2:')
x = amp(:,1); y = amp(:,2);
[H,P,CI,STATS] = ttest(x,y);
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,H,STATS] = signrank(x,y);
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);

disp('Ori1 vs Ori2 vs Plaid (Friedman''s test for K-related-samples):')
[P,TABLE,~] = friedman(amp(:,1:3),1,'off');
disp(['Chi^2(df=3' num2str(TABLE{2,3}) ') = ' num2str(TABLE{2,5},'%0.1f')]);
disp(['p            = ' num2str(P,'%0.3f')]);


% Compare delays
disp('-----------')
disp('Polar Delay')
disp('-----------')
disp('Ori vs Plaid (one-tail):')
x = delay(:,4); y = delay(:,3);
[H,P,CI,STATS] = ttest(x,y,'tail','right');
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,H,STATS] = signrank(x,y,'tail','right');
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,F] = circ_htest(x,y);
disp(' Hotelling''s test for angular means')
disp([' F=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P/2,'%0.2f')]);
disp('Ori1 vs Ori2:')
x = delay(:,1); y = delay(:,2);
[H,P,CI,STATS] = ttest(x,y);
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,H,STATS] = signrank(x,y);
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,F] = circ_htest(x,y);
disp(' Hotelling''s test for angular means')
disp([' F=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);

disp('Ori1 vs Plaid (one-tail):')
x = delay(:,1); y = delay(:,3);
[H,P,CI,STATS] = ttest(x,y,'tail','right');
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,H,STATS] = signrank(x,y,'tail','right');
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,F] = circ_htest(x,y);
disp(' Hotelling''s test for angular means')
disp([' F=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P/2,'%0.2f')]);

disp('Ori2 vs Plaid (one-tail):')
x = delay(:,2); y = delay(:,3);
[H,P,CI,STATS] = ttest(x,y,'tail','right');
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,H,STATS] = signrank(x,y,'tail','right');
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,F] = circ_htest(x,y);
disp(' Hotelling''s test for angular means')
disp([' F=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P/2,'%0.2f')]);

disp('Ori1 vs Ori2 vs Plais (Friedman''s test for K-related-samples):')
[P,TABLE,~] = friedman(delay(:,1:3),1,'off');
disp(['Chi^2(df=3' num2str(TABLE{2,3}) ') = ' num2str(TABLE{2,5},'%0.1f')]);
disp(['p            = ' num2str(P,'%0.3f')]);


%% Correlations
% Between-subject correlation
xData = [];
xLabel = [];
for subjInd = 1:length(subjList)
    for sessInd = 1:2

        tmp = squeeze(d{subjInd}.(['sess' num2str(sessInd)]).sin);
        tmp = mean(tmp,1);
        xData = cat(1,xData,tmp);
        sz = size(tmp);
        xLabel = cat(1,xLabel,[ones(sz(1),1).*subjInd ones(sz(1),1).*sessInd]);
    end
%     xData = cat(1,xData(1:end-2,:),mean(xData(end-1:end,:),1));
%     xLabel = cat(1,xLabel(1:end-2,:),xLabel(end,:));
end
xData = cat(2,mean(xData(:,1:2),2),xData(:,3));
xDataMean = mean(xData(:));

amp = abs(xData);
delay = angle(xData);
delay = wrapToPi(delay - angle(xDataMean));
delay = -(delay + angle(xDataMean))/pi*6;

fCorr1 = figure('WindowStyle','docked');
x = delay(:,2)-delay(:,1);
y = amp(:,2) - amp(:,1);
hScat = scatter(x,y);
[R,P] = corr(x,y,'type','Pearson','tail','left');
xlabel('delay (plaid - ori) in sec')
ylabel('amp (plaid - ori)')
refFit = polyfit(x,y,1);
title({'Subj*Sess' ['Pearson''s R=' num2str(R,'%0.2f') '; one-tailed p=' num2str(P,'%0.3f') '; slope=' num2str(refFit(1),'%0.2f')]})
refline
set(gca,'PlotBoxAspectRatio',[1 1 1])

xLabel = xLabel(:,1);
cMap = jet(length(subjList));
cData = nan(length(x),3);
for subjInd = 1:length(subjList)
    cData(xLabel==subjInd,:) = repmat(cMap(subjInd,:),[sum(xLabel==subjInd) 1]);
end
hScat.CData = cData;
hScat.MarkerFaceColor = hScat.MarkerEdgeColor;

if figOption.save
    writeFig(fCorr1,mfilename,'diffCorr_subjAndSess',verbose)
end


%average sessions
fCorr2 = figure('WindowStyle','docked');
delay = mean(cat(3,delay(1:2:end,:),delay(2:2:end,:)),3);
amp = mean(cat(3,amp(1:2:end,:),amp(2:2:end,:)),3);
x = delay(:,2)-delay(:,1);
y = amp(:,2) - amp(:,1);
hScat = scatter(x,y);
[R,P] = corr(x,y,'type','Spearman','tail','left');
xlabel('delay (plaid - ori) in sec')
ylabel('amp (plaid - ori)')
refFit = polyfit(x,y,1);
title({'Subj*Sess' ['Spearman''s R=' num2str(R,'%0.2f') '; one-tailed p=' num2str(P,'%0.3f') '; slope=' num2str(refFit(1),'%0.2f')]})
refline
set(gca,'PlotBoxAspectRatio',[1 1 1])
hScat.MarkerFaceColor = 'k';
hScat.MarkerEdgeColor = 'w';

if figOption.save
    writeFig(fCorr2,mfilename,'diffCorr_subj',verbose)
end

% Between-run correlation
%get grand mean
xData = [];
for subjInd = 1:length(subjList)
    for sessInd = 1:2
        for condInd = 1:3
            tmp = squeeze(d{subjInd}.(['sess' num2str(sessInd)]).sin);
            xData = [xData; tmp(:)];
        end
    end
end
xDataMean = mean(xData);
%remove subject*sess*cond effects
xData = [];
for subjInd = 1:length(subjList)
    for sessInd = 1:2
        for condInd = 1:3
            tmp = squeeze(d{subjInd}.(['sess' num2str(sessInd)]).sin);
            tmp = tmp - mean(tmp,1) + xDataMean;
            xData = [xData; tmp(:)];
        end
    end
end
amp = abs(xData);
delay = wrapToPi(angle(xData) - angle(xDataMean));
delay = -delay/pi*6;
delay = delay + (-angle(xDataMean)/pi*6);


fCorr3 = figure('WindowStyle','docked');
x = delay;
y = amp;
hScat = scatter(x,y);
[R,P] = corr(x,y,'tail','left');
xlabel('delay (sec)')
ylabel('amp (%BOLD)')
refFit = polyfit(x,y,1);
title({'Individual Runs, subj*sess*cond effect removed' ['Pearson''s R=' num2str(R,'%0.2f') '; one-tailed p=' num2str(P,'%0.3f') '; slope=' num2str(refFit(1),'%0.2f')]})
refline;
hScat.SizeData = hScat.SizeData/4;
hScat.MarkerFaceColor = 'k';
hScat.MarkerEdgeColor = 'w';
set(gca,'PlotBoxAspectRatio',[1 1 1])

if figOption.save
    writeFig(fCorr3,mfilename,'run2runCorr',verbose)
end

%% Reorder for notebook
fList = [fGroup([1 3 2]) fCorr1 fCorr2 fCorr3]';
delete(fList(4))
if ~verbose
    fList(4) = [];
end
set(0,'Children',flipud(fList))




