function exclusion = sinusoidalGroupAnalysis(threshType)
if ~exist('threshType','var') || isempty(threshType)
    threshType = 'p'; % 'none', 'p' or 'fdr'
end
doVein = 1;
noMovement = 1;
threshVal = 0.05;
plotAllSubj = 0;
saveFig = 0;

%colors
colors = [  0         0.4470    0.7410
            0.8500    0.3250    0.0980
            0.9290    0.6940    0.1250];
lw = 1.5;

if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
dataDir = 'C-derived\DecodingHR';
funPath = fullfile(repoPath,dataDir,'fun');
funLevel = 'z';
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';
% subjList = {'02jp' '03sk' '04sp'}';
if noMovement
    fileSuffix = '_maskSinAndHrFit_noMovement.mat';
else
    fileSuffix = '_maskSinAndHrFit.mat';
end

%make sure everything is forward slash for mac, linux pc compatibility
for tmpPath = {'repoPath' 'dataDir' 'funPath'}
    eval([char(tmpPath) '(strfind(' char(tmpPath) ',''\''))=''/'';']);
end



disp(['IN: Sinusoidal BOLD responses from anatomical V1 ROI (' fullfile(dataDir,funLevel) ')'])
disp(['threshVal=' num2str(threshVal)])
disp('F(IN)=OUT: threshold included voxels and analyse responses averaged across the ROI')
disp(['OUT: figures and stats'])



%% Load data
dAll = cell(size(subjList,1),1);
featSelStats = cell(size(subjList,1),1);
if doVein
    veinAll = cell(size(subjList,1),1);
end
for subjInd = 1:size(subjList,1)
    if doVein
        load(fullfile(funPath,funLevel,[subjList{subjInd} fileSuffix]),'d','featSelStats','vein');
        veinAll{subjInd} = vein;
    else
        load(fullfile(funPath,funLevel,[subjList{subjInd} fileSuffix]),'d','featSelStats');
    end
    dAll{subjInd} = d;
    featSelStatsAll{subjInd} = featSelStats;
end
d = dAll; clear dAll
featSelStats = featSelStatsAll; clear featSelStatsAll
if doVein
    vein = veinAll; clear veinAll
end

%% Process data
dP = d;

% Threshold and average voxels in cartesian space, using statistical tests
% on the other session
rLim = nan(length(dP),1);
for subjInd = 1:length(dP)
    if ~strcmp(threshType,'none')
        indAct = featSelStats{subjInd}.anyCondActivation.sess2.(upper(threshType))<threshVal;
        if doVein; indVein = vein{subjInd}.sess2.mask; end
        dP{subjInd}.sess1.data = dP{subjInd}.sess1.data(:,indAct&indVein,:,:);
        
        indAct = featSelStats{subjInd}.anyCondActivation.sess1.(upper(threshType))<threshVal;
        if doVein; indVein = vein{subjInd}.sess1.mask; end
        dP{subjInd}.sess2.data = dP{subjInd}.sess2.data(:,indAct&indVein,:,:);
    else
        if doVein; indVein = vein{subjInd}.sess2.mask; end
        dP{subjInd}.sess1.data = dP{subjInd}.sess1.data(:,indVein,:,:);
        if doVein; indVein = vein{subjInd}.sess1.mask; end
        dP{subjInd}.sess2.data = dP{subjInd}.sess2.data(:,indVein,:,:);
    end
    dP{subjInd}.sess1.data = mean(dP{subjInd}.sess1.data,2);
    dP{subjInd}.sess2.data = mean(dP{subjInd}.sess2.data,2);
    dP{subjInd}.sess1 = rmfield(dP{subjInd}.sess1,'hr');
    dP{subjInd}.sess2 = rmfield(dP{subjInd}.sess2,'hr');

    rLim(subjInd) = max(abs([dP{subjInd}.sess1.data(:); dP{subjInd}.sess2.data(:)]));
end

% Show scanner trigger problem
clear tmp tmp1 tmp2 tmpData
figure('WindowStyle','docked');
sz = 0;
for subjInd = 1:length(subjList)
    tmp = size(dP{subjInd}.sess1.data,1);
    if tmp>sz; sz = tmp; end
    tmpData(subjInd,1) = circ_mean(angle(dP{subjInd}.sess1.data(:)));
    tmpData(subjInd,2) = circ_mean(angle(dP{subjInd}.sess2.data(:)));
end

for subjInd = 1:length(subjList)
    tmp1 = squeeze(dP{subjInd}.sess1.data)';
    tmpInd = squeeze(dP{subjInd}.sess1.runLabel)';
    [~,b] = sort(tmpInd(:));
    tmpInd(b) = 1:numel(tmpInd);
    tmp1 = tmp1(tmpInd);
    tmp1 = tmp1(:);
    tmp1 = angle(tmp1);
    tmp1 = wrapToPi(tmp1-circ_mean(tmp1));
    
    tmp2 = squeeze(dP{subjInd}.sess2.data)';
    tmpInd = squeeze(dP{subjInd}.sess2.runLabel)';
    [~,b] = sort(tmpInd(:));
    tmpInd(b) = 1:numel(tmpInd);
    tmp2 = tmp2(tmpInd);
    tmp2 = tmp2(:);
    tmp2 = angle(tmp2);
    tmp2 = wrapToPi(tmp2-circ_mean(tmp2));
    
    tmp = cat(1,tmp1,nan(sz*3-length(tmp1)+1,1),tmp2);
%     tmp = wrapToPi(tmp-circ_mean(tmp(~isnan(tmp))));
    
%     tmp = wrapToPi(tmp + circ_mean(tmpData(:)));
    
    hp(subjInd) = plot(tmp,'-o'); hold on
end
TR = 1;
stimCycleLength = 12;
TRrad = TR/stimCycleLength*(2*pi);
ax = gca;
yTickRad = -pi:TRrad:pi;
yTickSec = yTickRad/(2*pi)*stimCycleLength;
ax.YTick = yTickRad;
ax.YTickLabel = num2str(yTickSec');
ylim([yTickRad(1) yTickRad(end)])
ax.YGrid = 'on';
xlabel({'Functional runs' 'in order of acquisition'})
ylabel('TRs');
legend(char(subjList))

% Exclude
subjInd = 2;
sessInd = 1;
dataTmp = angle(dP{subjInd}.(['sess' num2str(sessInd)]).data);
dataTmp = wrapToPi(dataTmp - circ_mean(dataTmp));
dataTmp = abs(dataTmp);
a = max(dataTmp(:),[],1);
[runInd,condInd] = find(squeeze(dataTmp)==a);
exclusion.subj = subjInd;
exclusion.sess = sessInd;
exclusion.run = runInd;
exclusion.cond = condInd;
% exclusion.subj = subjList(subjInd)
% exclusion.subj = {'03sk'};
% exclusion.sess = 1;
% exclusion.run = 4;
% exclusion.cond = 2;

if exist('exclusion','var') && ~isempty(exclusion) && ~isempty(exclusion.subj)
    for i = 1:length(exclusion.subj)
        dP{exclusion.subj(i)}.(['sess' num2str(i)]).data(exclusion.run,:,:) = [];
        dP{exclusion.subj(i)}.(['sess' num2str(i)]).runLabel(exclusion.run,:,:) = [];
    end
end


% Plot single subjects
for subjInd = 1:length(subjList)
    if plotAllSubj || subjInd==1
        fSubj(subjInd) = figure('WindowStyle','docked');
        for sessInd = 1:2
            subplot(1,2,sessInd)
            xData = squeeze(dP{subjInd}.(['sess' num2str(sessInd)]).data);
            for condInd = 1:size(xData,2)
                [theta,rho] = cart2pol(real(xData(:,condInd)),imag(xData(:,condInd)));
                h(condInd) = polarplot(theta,rho,'o'); hold on
                h(condInd).Color = colors(condInd,:);
                [theta,rho] = cart2pol(real(mean(xData(:,condInd),1)),imag(mean(xData(:,condInd),1)));
                hAv(condInd) = polarplot(theta,rho,'o'); hold on
                hAv(condInd).MarkerEdgeColor = 'w';
                hAv(condInd).MarkerFaceColor = colors(condInd,:);
                hAv(condInd).LineWidth = 0.25;
            end
            uistack(hAv,'top')
            title(['Sess' num2str(sessInd)])
            rlim([0 max(rLim)])
            
            ax = gca;
            ax.ThetaTickLabel = 12-ax.ThetaTick(1:end)/360*12;
            ax.ThetaTickLabel(1,:) = '0 ';
%             ax.ThetaAxis.Label.String = {'delay' '(sec)'};
%             ax.ThetaAxis.Label.Rotation = 0;
%             ax.ThetaAxis.Label.HorizontalAlignment = 'left';
%             ax.RAxis.Label.String = 'amp (%BOLD)';
%             ax.RAxis.Label.Rotation = 80;
        end
        suptitle(subjList{subjInd})
        hl = legend(char({'ori1' 'ori2' 'plaid'}),'Location','east');
        hl.Color = 'none';
        hl.Box = 'off';
        
        if saveFig
            filename = fullfile(pwd,mfilename);
            if ~exist(filename,'dir'); mkdir(filename); end
            filename = fullfile(filename,subjList{subjInd});
            fSubj(subjInd).Color = 'none';
            set(findobj(fSubj(subjInd).Children,'type','Axes'),'color','none')
            set(findobj(fSubj(subjInd).Children,'type','PolarAxes'),'color','none')
            saveas(fSubj(subjInd),[filename '.svg']); disp([filename '.svg'])
            fSubj(subjInd).Color = 'w';
%             set(findobj(fSubj(subjInd).Children,'type','Axes'),'color','w')
%             set(findobj(fSubj(subjInd).Children,'type','PolarAxes'),'color','w')
            saveas(fSubj(subjInd),filename); disp([filename '.fig'])
            saveas(fSubj(subjInd),filename); disp([filename '.jpg'])
        end
    end
end


% Average runs in cartesian space
condList = {'ori1' 'ori2' 'plaid'};
xData = nan(length(subjList),length(condList),2);
for subjInd = 1:length(subjList)
    xData(subjInd,:,:) = permute(cat(1,mean(dP{subjInd}.sess1.data,1),mean(dP{subjInd}.sess2.data,1)),[2 3 1]);
end
xDataInfo = 'subj x cond[or1, ori2, plaid] x sess';

% Average ori cartesian space
xData = cat(2,xData,mean(xData(:,1:2,:),2));
xDataInfo = 'subj x cond[or1, ori2, plaid, ori] x sess';

% Average sessions in cartesian space
xData = mean(xData,3);
xDataInfo = 'subj x cond[or1, ori2, plaid, ori]';

% Save intermediary data for hrGroupAnalysis
save(fullfile(funPath,funLevel,'tmp.mat'),'xData','xDataInfo')

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

if saveFig
    filename = fullfile(pwd,mfilename);
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,'group');
    fGroup(1).Color = 'none';
    set(findobj(fGroup(1).Children,'type','Axes'),'color','none')
    set(findobj(fGroup(1).Children,'type','PolarAxes'),'color','none')
    saveas(fGroup(1),[filename '.svg']); disp([filename '.svg'])
    fGroup(1).Color = 'w';
%     set(findobj(fGroup(1).Children,'type','Axes'),'color','w')
%     set(findobj(fGroup(1).Children,'type','PolarAxes'),'color','w')
    saveas(fGroup(1),filename); disp([filename '.fig'])
    saveas(fGroup(1),filename); disp([filename '.jpg'])
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

if saveFig
    filename = fullfile(pwd,mfilename);
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,'groupAmp');
    fGroup(2).Color = 'none';
    set(findobj(fGroup(2).Children,'type','Axes'),'color','none')
    saveas(fGroup(2),[filename '.svg']); disp([filename '.svg'])
    fGroup(2).Color = 'w';
%     set(findobj(fGroup(2).Children,'type','Axes'),'color','w')
    saveas(fGroup(2),filename); disp([filename '.fig'])
    saveas(fGroup(2),filename); disp([filename '.jpg'])
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

disp('***')
disp(['delay (plaid-ori) = ' num2str(diff(xData_thetaMean(:,[4 3])),'%0.3fs')])
disp(['delay (ori2-ori1) = ' num2str(diff(xData_thetaMean(:,[1 2])),'%0.3fs')])
disp(['amp (plaid-ori) = ' num2str(diff(xData_rhoMean(:,[4 3])),'%0.3f%%BOLD')])
disp(['amp (ori2-ori1) = ' num2str(diff(xData_rhoMean(:,[1 2])),'%0.3f%%BOLD')])
disp('***')

if saveFig
    filename = fullfile(pwd,mfilename);
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,'groupDelay');
    fGroup(3).Color = 'none';
    set(findobj(fGroup(3).Children,'type','Axes'),'color','none')
    saveas(fGroup(3),[filename '.svg']); disp([filename '.svg'])
    fGroup(3).Color = 'w';
%     set(findobj(fGroup(3).Children,'type','Axes'),'color','w')
    saveas(fGroup(3),filename); disp([filename '.fig'])
    saveas(fGroup(3),filename); disp([filename '.jpg'])
end

%% Stats
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
x = abs(xData(:,4)); y = abs(xData(:,3));
[H,P,CI,STATS] = ttest(x,y,'tail','right');
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,H,STATS] = signrank(x,y,'tail','right');
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);
disp('Ori1 vs Ori2:')
x = abs(xData(:,1)); y = abs(xData(:,2));
[H,P,CI,STATS] = ttest(x,y);
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,H,STATS] = signrank(x,y);
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);

% Compare delays
disp('-----------')
disp('Polar Delay')
disp('-----------')
disp('Ori vs Plaid (one-tail):')
x = angle(xData(:,4)); y = angle(xData(:,3));
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
x = angle(xData(:,1)); y = angle(xData(:,2));
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
x = angle(xData(:,1)); y = angle(xData(:,3));
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
x = angle(xData(:,2)); y = angle(xData(:,3));
[H,P,CI,STATS] = ttest(x,y,'tail','right');
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,H,STATS] = signrank(x,y,'tail','right');
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,F] = circ_htest(x,y);
disp(' Hotelling''s test for angular means')
disp([' F=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P/2,'%0.2f')]);


%% Correlations
% Between-subject correlation
xData = [];
xLabel = [];
for subjInd = 1:length(subjList)
    for sessInd = 1:2
        tmp = squeeze(dP{subjInd}.(['sess' num2str(sessInd)]).data);
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

f = figure('WindowStyle','docked');
x = wrapToPi(delay(:,2)-delay(:,1));
y = amp(:,2) - amp(:,1);
hScat = scatter(x,y);
[R,P] = corr(x,y,'type','Spearman');
xlabel('delay (plaid - ori)')
ylabel('amp (plaid - ori)')
refFit = polyfit(x,y,1);
title({'Subj*Sess' ['Spearman''s R=' num2str(R,'%0.2f') '; p=' num2str(P,'%0.3f') '; slope=' num2str(refFit(1),'%0.2f')]})
refline


xLabel = xLabel(:,1);
cMap = jet(length(subjList));
cData = nan(length(x),3);
for subjInd = 1:length(subjList)
    cData(xLabel==subjInd,:) = repmat(cMap(subjInd,:),[sum(xLabel==subjInd) 1]);
end
hScat.CData = cData;
hScat.MarkerFaceColor = hScat.MarkerEdgeColor;


% Between-run correlation
xData = [];
xLabel = [];
i = 0;
for subjInd = 1:length(subjList)
    for sessInd = 1:2
        i = i+1;
        tmp = squeeze(dP{subjInd}.(['sess' num2str(sessInd)]).data);
        xData = cat(1,xData,tmp);
        sz = size(tmp);
        xLabel = cat(1,xLabel,[ones(sz(1),1).*subjInd ones(sz(1),1).*sessInd ones(sz(1),1).*i]);
    end
end
% remove subject*cond effects
xDataMean = mean(xData(:));
for ii = 1:length(unique(xLabel(:,end)))
    ind = ii==xLabel(:,end);
    xData(ind,:) = xData(ind,:) - mean(xData(ind,:),1);
end
xData = xData+xDataMean;

amp = abs(xData);
delay = angle(xData);
delay = wrapToPi(delay - angle(xDataMean));
delay = -(delay + angle(xDataMean))/pi*6;

f = figure('WindowStyle','docked');
x = delay(:);
y = amp(:);
hScat = scatter(x,y);
[R,P] = corr(x,y);
xlabel('delay')
ylabel('amp')
refFit = polyfit(x,y,1);
title({'Individual Runs, subj*sess*cond effect removed' ['Pearson''s R=' num2str(R,'%0.2f') '; p=' num2str(P,'%0.3f') '; slope=' num2str(refFit(1),'%0.2f')]})
refline;

