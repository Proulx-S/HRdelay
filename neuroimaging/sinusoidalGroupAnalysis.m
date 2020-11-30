function sinusoidalGroupAnalysis(threshType)
if ~exist('threshType','var')
    threshType = 'fdr'; % 'none', 'p' or 'fdr'
end
threshVal = 0.05;
adjVoxDelay = 0;

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
funLevel = 'zSin';
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';
fileSuffix = '_maskSinAndHrFit.mat';

%make sure everything is forward slash for mac, linux pc compatibility
for tmpPath = {'repoPath' 'dataDir' 'funPath'}
    eval([char(tmpPath) '(strfind(' char(tmpPath) ',''\''))=''/'';']);
end



disp(['IN: Sinusoidal BOLD responses from anatomical V1 ROI (' fullfile(dataDir,funLevel) ')'])
disp(['threshVal=' num2str(threshVal)])
disp(['adjVoxDelay=' num2str(adjVoxDelay)])
disp('F(IN)=OUT: threshold included voxels and analyse responses averaged across the ROI')
disp(['OUT: figures and stats'])



%% Load data
dAll = cell(size(subjList,1),1);
for subjInd = 1:size(subjList,1)
    load(fullfile(funPath,funLevel,[subjList{subjInd} fileSuffix]),'d');
    dAll{subjInd} = d;
end
d = dAll; clear dAll

%% Process data
dP = d;
% Remove inter-voxel variations in delay, using delay from the other
% session
if adjVoxDelay
    voxDelaySess1 = angle(mean(mean(dP{subjInd}.sess1.xData,3),1)) - angle(mean(dP{subjInd}.sess1.xData(:)));
    voxDelaySess2 = angle(mean(mean(dP{subjInd}.sess2.xData,3),1)) - angle(mean(dP{subjInd}.sess2.xData(:)));
    
    theta = angle(dP{subjInd}.sess1.xData) - voxDelaySess2;
    rho = abs(dP{subjInd}.sess1.xData);
    [X,Y] = pol2cart(theta,rho);
    dP{subjInd}.sess1.xData = complex(X,Y);
    
    theta = angle(dP{subjInd}.sess2.xData) - voxDelaySess1;
    rho = abs(dP{subjInd}.sess2.xData);
    [X,Y] = pol2cart(theta,rho);
    dP{subjInd}.sess2.xData = complex(X,Y); clear X Y theta rho voxDelaySess1 voxDelaySess2
end

% figure('WindowStyle','docked')
% subplot(2,1,1)
% % polarplot(squeeze(mean(d{subjInd}.sess1.xData,1)),'o');
% polarplot(squeeze(mean(d{subjInd}.sess1.xData,2)),'o');
% abs(mean(d{subjInd}.sess1.xData(:)))
% subplot(2,1,2)
% % polarplot(squeeze(mean(dP{subjInd}.sess1.xData,1)),'o');
% polarplot(squeeze(mean(dP{subjInd}.sess1.xData,2)),'o');
% abs(mean(dP{subjInd}.sess1.xData(:)))

% Threshold and average voxels in cartesian space, using statistical tests
% on the other session
rLim = (1:length(dP))';
for subjInd = 1:length(dP)
    switch threshType
        % no voxel selection
        case 'none' 
            dP{subjInd}.sess1.xData = mean(dP{subjInd}.sess1.xData,2);
            dP{subjInd}.sess2.xData = mean(dP{subjInd}.sess2.xData,2);
        % voxel selection cross-validated between sessions
        case 'p'
            dP{subjInd}.sess1.xData = mean(dP{subjInd}.sess1.xData(:,dP{subjInd}.sess2.P<threshVal,:),2);
            dP{subjInd}.sess2.xData = mean(dP{subjInd}.sess2.xData(:,dP{subjInd}.sess1.P<threshVal,:),2);
        case 'fdr'
            dP{subjInd}.sess1.xData = mean(dP{subjInd}.sess1.xData(:,dP{subjInd}.sess2.FDR<threshVal,:),2);
            dP{subjInd}.sess2.xData = mean(dP{subjInd}.sess2.xData(:,dP{subjInd}.sess1.FDR<threshVal,:),2);
        otherwise
            error('X')
    end
    
    rLim(subjInd) = max(abs([dP{subjInd}.sess1.xData(:); dP{subjInd}.sess2.xData(:)]));
end

% Plot single subjects
plotAll = 0;
for subjInd = 1:length(subjList)
    if plotAll || subjInd==1
        fSubj(subjInd) = figure('WindowStyle','docked');
        for sessInd = 1:2
            subplot(1,2,sessInd)
            xData = squeeze(dP{subjInd}.(['sess' num2str(sessInd)]).xData);
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
        
        filename = fullfile(pwd,mfilename);
        if ~exist(filename,'dir'); mkdir(filename); end
        filename = fullfile(filename,subjList{subjInd});
        fSubj(subjInd).Color = 'none';
        set(findobj(fSubj(subjInd).Children,'type','Axes'),'color','none')
        set(findobj(fSubj(subjInd).Children,'type','PolarAxes'),'color','none')
        saveas(fSubj(subjInd),[filename '.svg']); disp([filename '.svg'])
        fSubj(subjInd).Color = 'w';
        set(findobj(fSubj(subjInd).Children,'type','Axes'),'color','w')
        set(findobj(fSubj(subjInd).Children,'type','PolarAxes'),'color','w')
        saveas(fSubj(subjInd),filename); disp([filename '.fig'])
    end
end

% Average runs in cartesian space
condList = {'ori1' 'ori2' 'plaid'};
xData = nan(length(subjList),length(condList),2);
for subjInd = 1:length(subjList)
    xData(subjInd,:,:) = permute(cat(1,mean(dP{subjInd}.sess1.xData,1),mean(dP{subjInd}.sess2.xData,1)),[2 3 1]);
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

filename = fullfile(pwd,mfilename);
if ~exist(filename,'dir'); mkdir(filename); end
filename = fullfile(filename,'group');
fGroup(1).Color = 'none';
set(findobj(fGroup(1).Children,'type','Axes'),'color','none')
set(findobj(fGroup(1).Children,'type','PolarAxes'),'color','none')
saveas(fGroup(1),[filename '.svg']); disp([filename '.svg'])
fGroup(1).Color = 'w';
set(findobj(fGroup(1).Children,'type','Axes'),'color','w')
set(findobj(fGroup(1).Children,'type','PolarAxes'),'color','w')
saveas(fGroup(1),filename); disp([filename '.fig'])


% Polar space
xData_theta = -angle(xData)./pi*6;
xData_rho = abs(xData);
xDataMean = mean(xData,1);
xDataMean_theta = -angle(xDataMean)./pi*6;
xDataMean_rho = abs(xDataMean);

fGroup(2) = figure('WindowStyle','docked');
[ax, pos] = tight_subplot(1,2);
spInd = 1;
axes(ax(spInd));
hb = bar(xDataMean_rho(:,[4 3])); hold on
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
hb = bar(xDataMean_rho(:,[1 2])); hold on
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
disp('***')
disp(['delay diff=' num2str(abs(diff(xDataMean_rho(:,[1 2]))),'%0.3fms')])
disp('***')

filename = fullfile(pwd,mfilename);
if ~exist(filename,'dir'); mkdir(filename); end
filename = fullfile(filename,'groupAmp');
fGroup(2).Color = 'none';
set(findobj(fGroup(2).Children,'type','Axes'),'color','none')
saveas(fGroup(2),[filename '.svg']); disp([filename '.svg'])
fGroup(2).Color = 'w';
set(findobj(fGroup(2).Children,'type','Axes'),'color','w')
saveas(fGroup(2),filename); disp([filename '.fig'])


fGroup(3) = figure('WindowStyle','docked');
[ax, pos] = tight_subplot(2, 1,[0],[0.1 0],[0.1 0]);
spInd = 1;
axes(ax(spInd));
hb = barh(xDataMean_theta(:,[4 3])); hold on
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
hb = barh(xDataMean_theta(:,[1 2])); hold on
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

filename = fullfile(pwd,mfilename);
if ~exist(filename,'dir'); mkdir(filename); end
filename = fullfile(filename,'groupDelay');
fGroup(3).Color = 'none';
set(findobj(fGroup(3).Children,'type','Axes'),'color','none')
saveas(fGroup(3),[filename '.svg']); disp([filename '.svg'])
fGroup(3).Color = 'w';
set(findobj(fGroup(3).Children,'type','Axes'),'color','w')
saveas(fGroup(3),filename); disp([filename '.fig'])


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
disp('Ori vs Plaid:')
x = abs(xData(:,4)); y = abs(xData(:,3));
[H,P,CI,STATS] = ttest(x,y);
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,H,STATS] = signrank(x,y);
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
disp('Ori vs Plaid:')
x = angle(xData(:,4)); y = angle(xData(:,3));
[H,P,CI,STATS] = ttest(x,y);
disp(' Student''s t-test')
disp([' t=' num2str(STATS.tstat,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,H,STATS] = signrank(x,y);
disp(' Wilcoxon signed rank test')
disp([' signed rank=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);
[P,F] = circ_htest(x,y);
disp(' Hotelling''s test for angular means')
disp([' F=' num2str(STATS.signedrank,'%0.2f') '; p=' num2str(P,'%0.2f')]);
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

