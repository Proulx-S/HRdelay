clear all
close all


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
xData_theta = angle(xData);
xData_rho = abs(xData);

% Plot
xDataNorm = xData - mean(xData(:,1:3),2) + mean(mean(xData(:,1:3),2),1);
xDataNorm_theta = angle(xData);
xDataNorm_rho = abs(xData);
figure('WindowStyle','docked');
h = polarplot(xDataNorm_theta(:,1:3),xDataNorm_rho(:,1:3),'o'); hold on
bar(xData_rho)

% Compare amplitudes
plaidInd = 3;
condInd = 4; % 1; 2; [1 2];
disp('amplitude')
[H,P,CI,STATS] = ttest(mean(xData_rho(condInd,:),1),xData_rho(plaidInd,:));
[STATS.tstat P]
[P,H,STATS] = signrank(mean(xData_rho(condInd,:),1),xData_rho(plaidInd,:));
[STATS.signedrank P]
% Compare delays
disp('delay')
[H,P,CI,STATS] = ttest(mean(xData_theta(condInd,:),1),xData_theta(plaidInd,:));
[STATS.tstat P]
[P,H,STATS] = signrank(mean(xData_theta(condInd,:),1),xData_theta(plaidInd,:));
[STATS.signedrank P]
[P,F] = circ_htest(mean(xData_theta(condInd,:),1),xData_theta(plaidInd,:));
[F P]




















load(fullfile(repoPath,dataDir,funPath,funLevel,[subjList{subjInd} '.mat']),'d')
ori1{subjInd,sessInd} = d.xData(d.label==1,:);
ori2{subjInd,sessInd} = d.xData(d.label==2,:);
plaid{subjInd,sessInd} = d.normData;


data.ori1 = ori1;
data.ori2 = ori2;
data.plaid = plaid;


data = loadData(dataRepo,dataDir,dataLevel,fileList)






dataRepo_in = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
dataRepo_out = dataRepo_in;
dataDir = 'C-derived\DecodingHR';
dataLevel = 'z';
fileList = {'02jp_sess1' '03sk_sess1' '04sp_sess1' '05bm_sess1' '06sb_sess1' '07bj_sess1';...
    '02jp_sess2' '03sk_sess2' '04sp_sess2' '05bm_sess2' '06sb_sess2' '07bj_sess2'}';
fileList


% Load complex numbers representing the sinusoidal response
data = loadData(dataRepo,dataDir,dataLevel,fileList)
data.ori1

% Average voxels in cartesian space
dataP = avVox(data);
dataP.ori1
rLim = getRlim(dataP);

% Plot single subjects
plotSubj(dataP,1,1);
% for subjInd = 1:6
%     plotSubj(dataP,subjInd,1,rLim);
%     plotSubj(dataP,subjInd,2,rLim);
% end

% % Remove session effects (in cartesian space)
% sessMeans = mean(getSessMeans(dataP),1);
% subjMeans = mean(sessMeans,3);
% condList = fields(dataP);
% for condInd = 1:length(condList)
%     for subjInd = 1:6
%         for sessInd = 1:2
%             dataP.(condList{condInd}){subjInd,sessInd} = dataP.(condList{condInd}){subjInd,sessInd} - sessMeans(:,subjInd,sessInd);
%             dataP.(condList{condInd}){subjInd,sessInd} = dataP.(condList{condInd}){subjInd,sessInd} + subjMeans(:,subjInd);
%         end
%     end
% end

% % Plot single subjects
% for subjInd = 1:6
%     plotSubj(dataP,subjInd,1,rLim);
%     plotSubj(dataP,subjInd,2,rLim);
% end

% Get session means (in cartesian space)
sessMeans = getSessMeans(dataP); % cond[ori1 ori2 plaid ori] x subj x sess

% Get subject means in cartesian space then move to polar space
subjMeans = mean(sessMeans,3);

% Move to polar space
xData_theta = angle(subjMeans);
xData_rho = abs(subjMeans);


% Compare amplitudes
condInd = 4; % 1; 2; [1 2];
disp('amplitude')
[H,P,CI,STATS] = ttest(mean(xData_rho(condInd,:),1),xData_rho(3,:));
[STATS.tstat P]
[P,H,STATS] = signrank(mean(xData_rho(condInd,:),1),xData_rho(3,:));
[STATS.signedrank P]
% Compare delays
disp('delay')
[H,P,CI,STATS] = ttest(mean(xData_theta(condInd,:),1),xData_theta(3,:));
[STATS.tstat P]
[P,H,STATS] = signrank(mean(xData_theta(condInd,:),1),xData_theta(3,:));
[STATS.signedrank P]
[P,F] = circ_htest(mean(xData_theta(condInd,:),1),xData_theta(3,:));
[F P]

%
%
% % ori1
% %% Compare amplitudes
% disp('amplitude')
% [H,P,CI,STATS] = ttest(abs(subjMeans(1,:)),abs(subjMeans(3,:)));
% P
% [P,H,STATS] = signrank(abs(subjMeans(1,:)),abs(subjMeans(3,:)));
% P
% %% Compare delays
% disp('delay')
% [H,P,CI,STATS] = ttest(angle(subjMeans(1,:)),angle(subjMeans(3,:)));
% P
% [P,H,STATS] = signrank(angle(subjMeans(1,:)),angle(subjMeans(3,:)));
% P
% [P,F] = circ_htest(angle(subjMeans(1,:)),angle(subjMeans(3,:)));
% P
