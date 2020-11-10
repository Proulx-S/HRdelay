function sinusoidalGroupAnalysis

repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
dataDir = 'C-derived\DecodingHR';
funPath = fullfile(repoPath,dataDir,'fun');
funLevel = 'zSin';
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';
fileSuffix = '_maskSinFit.mat';



%% Load data
dAll = cell(size(subjList,1),1);
for subjInd = 1:size(subjList,1)
    load(fullfile(funPath,funLevel,[subjList{subjInd} fileSuffix]),'d');
    dAll{subjInd} = d;
end
d = dAll; clear dAll

%% Process data
% Average voxels in cartesian space
dP = d;
rLim = (1:length(d))';
for subjInd = 1:length(d)
    dP{subjInd}.sess1.xData = mean(dP{subjInd}.sess1.xData,2);
    dP{subjInd}.sess2.xData = mean(dP{subjInd}.sess2.xData,2);
    
    rLim(subjInd) = max(abs([dP{subjInd}.sess1.xData(:); dP{subjInd}.sess2.xData(:)]));
end

% Plot single subjects
plotAll = 1;
for subjInd = 1:length(subjList)
    if plotAll || subjInd==1
        figure('WindowStyle','docked');
        for sessInd = 1:2
            subplot(1,2,sessInd)
            xData = squeeze(dP{subjInd}.(['sess' num2str(sessInd)]).xData);
            for condInd = 1:size(xData,2)
                [theta,rho] = cart2pol(real(xData(:,condInd)),imag(xData(:,condInd)));
                h(condInd) = polarplot(theta,rho,'o'); hold on
                [theta,rho] = cart2pol(real(mean(xData(:,condInd),1)),imag(mean(xData(:,condInd),1)));
                hAv(condInd) = polarplot(theta,rho,'o'); hold on
                hAv(condInd).MarkerEdgeColor = 'w';
                hAv(condInd).MarkerFaceColor = h(condInd).Color;
                hAv(condInd).LineWidth = 0.25;
            end
            uistack(hAv,'top')
            title(['Sess' num2str(sessInd)])
            rlim([0 max(rLim)])
        end
        suptitle(subjList{subjInd})
        hl = legend(char({'ori1' 'ori2' 'plaid'}),'Location','east');
        hl.Color = 'none';
        hl.Box = 'off';
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

% Move to polar space

%% Plot results
% Cartesian space
xDataNorm = xData - mean(xData(:,1:3),2) + mean(mean(xData(:,1:3),2),1);
xDataNormMean = mean(xDataNorm,1);
xDataNorm_theta = angle(xDataNorm);
xDataNorm_rho = abs(xDataNorm);
xDataNormMean_theta = angle(xDataNormMean);
xDataNormMean_rho = abs(xDataNormMean);
figure('WindowStyle','docked');
for condInd = 1:3
    h(condInd) = polarplot(xDataNorm_theta(:,condInd),xDataNorm_rho(:,condInd),'o'); hold on
    hM(condInd) = polarplot(xDataNormMean_theta(:,condInd),xDataNormMean_rho(:,condInd),'o'); hold on
    hM(condInd).MarkerFaceColor = h(condInd).Color;
    hM(condInd).MarkerEdgeColor = 'w';
end
uistack(hM,'top')
hl = legend(char({'ori1' 'ori2' 'plaid'}),'Location','east');
hl.Color = 'none';
hl.Box = 'off';
title('Group Response')

% Polar space
xData_theta = angle(xData);
xData_rho = abs(xData);
xDataMean = mean(xData,1);
xDataMean_theta = angle(xDataMean);
xDataMean_rho = abs(xDataMean);

figure('WindowStyle','docked');
hb = bar(xDataMean_rho(:,[4 3])); hold on
hp = plot(hb.XData,xData_rho(:,[4 3])','k');
xlim([0.25 2.75])
ylabel({'Response amplitude' '(%BOLD)'})
ax = gca; ax.XTickLabel = {'ori' 'plaid'};
ylim([0 max(xData_rho(:)).*1.1])
ax.PlotBoxAspectRatio = [0.2 1 1];
figure('WindowStyle','docked');
hb = bar(xDataMean_rho(:,[1 2])); hold on
hp = plot(hb.XData,xData_rho(:,[1 2])','k');
xlim([0.25 2.75])
ylabel({'Response amplitude' '(%BOLD)'})
ax = gca; ax.XTickLabel = {'ori1' 'ori2'};
ylim([0 max(xData_rho(:)).*1.1])
ax.PlotBoxAspectRatio = [0.2 1 1];

figure('WindowStyle','docked');
hb = barh(xDataMean_theta(:,[4 3])); hold on
hp = plot(xData_theta(:,[4 3])',hb.XData,'k');
ylim([0.25 2.75]);
xlabel({'Response Delay' '(rad)'})
ax = gca; ax.YTickLabel = {'ori' 'plaid'};
xlim([min(xData_theta(:)).*1.1 0])
ax.PlotBoxAspectRatio = [1 0.1 1];
figure('WindowStyle','docked');
hb = barh(xDataMean_theta(:,[1 2])); hold on
hp = plot(xData_theta(:,[1 2])',hb.XData,'k');
ylim([0.25 2.75]);
xlabel({'Response Delay' '(rad)'})
ax = gca; ax.YTickLabel = {'ori1' 'ori2'};
xlim([min(xData_theta(:)).*1.1 0])
ax.PlotBoxAspectRatio = [1 0.1 1];

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

