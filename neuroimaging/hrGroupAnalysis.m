function hrGroupAnalysis(threshType)
if ~exist('threshType','var')
    threshType = 'none'; % 'none', 'p' or 'fdr'
end
threshVal = 0.05;
plotAllSubj = 0;


repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
dataDir = 'C-derived\DecodingHR';
funPath = fullfile(repoPath,dataDir,'fun');
funLevel = 'zHr';
funLevelSin = 'zSin';
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';
fileSuffix = '_maskSinAndHrFit.mat';

disp(['IN: BOLD hemodynamic responses (HR) from anatomical V1 ROI (' fullfile(dataDir,funLevel) ')'])
disp('F(IN)->OUT: threshold included voxels and analyse HR averaged across the ROI')
disp(['OUT: figures'])



%% Load data
hrAll = cell(size(subjList,1),1);
if ~strcmp(threshType,'none')
    dForThreshAll = cell(size(subjList,1),1);
end
for subjInd = 1:size(subjList,1)
    load(fullfile(funPath,funLevel,[subjList{subjInd} fileSuffix]),'hr');
    hrAll{subjInd} = hr;
    switch threshType
        case 'none'
        case 'p'
            load(fullfile(funPath,'zSin',[subjList{subjInd} fileSuffix]),'d');
            dForThreshAll{subjInd}.sess1 = d.sess1.P;
            dForThreshAll{subjInd}.sess2 = d.sess2.P;
        case 'fdr'
            load(fullfile(funPath,'zSin',[subjList{subjInd} fileSuffix]),'d');
            dForThreshAll{subjInd}.sess1 = d.sess1.FDR;
            dForThreshAll{subjInd}.sess2 = d.sess2.FDR;
        otherwise
            error('X')
    end
end
hr = hrAll; clear dAll
if ~strcmp(threshType,'none')
    dForThresh = dForThreshAll; clear dForThreshAll
end
%% Process data
% Average voxels
hrP = hr;
rLim = (1:length(hr))';
for subjInd = 1:length(hr)
    switch threshType
        % no voxel selection
        case 'none' 
            hrP{subjInd}.sess1 = mean(hrP{subjInd}.sess1,2);
            hrP{subjInd}.sess2 = mean(hrP{subjInd}.sess2,2);
        % voxel selection cross-validated between sessions
        case {'p' 'fdr'}
            hrP{subjInd}.sess1 = mean(hrP{subjInd}.sess1(:,dForThresh{subjInd}.sess2<threshVal,:,:),2);
            hrP{subjInd}.sess2 = mean(hrP{subjInd}.sess2(:,dForThresh{subjInd}.sess1<threshVal,:,:),2);
        otherwise
            error('X')
    end
    rLim(subjInd) = max(abs([hrP{subjInd}.sess1(:); hrP{subjInd}.sess2(:)]));
end

% Plot single subjects
for subjInd = 1:length(subjList)
    if plotAllSubj || subjInd==1
        fSubj(subjInd) = figure('WindowStyle','docked');
        for sessInd = 1:2
            subplot(1,2,sessInd)
            xData = squeeze(hrP{subjInd}.(['sess' num2str(sessInd)]));
            for condInd = 1:size(xData,2)
                y = squeeze(xData(:,condInd,:));
                t = 0:11;
                yMean = mean(y,1);
                ySEM = std(y,[],1)./sqrt(size(y,1));
                he = errorbar(t,yMean,ySEM); hold on
                he.CapSize = 0;
            end
            title(['Sess' num2str(sessInd)])
            ylim([min(xData(:)) max(xData(:))])
        end
        legend(char({'ori1 +/-SEM' 'ori2 +/-SEM' 'plaid +/-SEM'}),'Location','south','box','off','color','none')
        suptitle(subjList{subjInd})
    end
end
clear xData

% Remove subjet effect
sinData = load(fullfile(funPath,funLevelSin,'tmp.mat'),'xData','xDataInfo');
rhoGroup = abs(mean(mean(sinData.xData(:,1:3),2),1));
thetaGroup = angle(mean(mean(sinData.xData(:,1:3),2),1));
t = 0:11;
for subjInd = 1:length(subjList)
    rhoSubj = abs(mean(sinData.xData(subjInd,1:3),2));
    thetaSubj = angle(mean(sinData.xData(subjInd,1:3),2));
    
    % response amplitude
    hrP{subjInd}.sess1 = hrP{subjInd}.sess1./rhoSubj.*rhoGroup;
    hrP{subjInd}.sess2 = hrP{subjInd}.sess2./rhoSubj.*rhoGroup;
    
    % response delay
    interpMethod = 'cubic'; % cubic convolution as implemented in Matlab2020b
    for sessInd = 1:2
        %upsample
        curHR = hrP{subjInd}.(['sess' num2str(sessInd)]);
        curT = permute(linspace(0,11,size(curHR,4)),[1 3 4 2]);

        n = size(curHR,4);
        curHR = permute(curHR,[4 1 2 3]);
        curT = permute(curT,[4 1 2 3]);
        curHR2 = cat(1,curHR              ,curHR              ,curHR              ,curHR              ,curHR);
        curT2 =  cat(1,curT-2*n,curT-1*n,curT-0*n,curT+1*n,curT+2*n);
        
        delta = 0.001;
        curT3 = (min(curT2):delta:(max(curT2)+1-delta))';
        curHR3 = nan([length(curT3) size(curHR2,[2 3 4])]);
        for i = 1:prod(size(curHR2,[2 3 4]))
            curHR3(:,i) = interp1(curT2,curHR2(:,i),curT3,interpMethod);
        end
%         figure('WindowStyle','docked');
%         plot(curT2,curHR2(:,i),'o'); hold on
%         plot(curT3,curHR3(:,i),'-');
        
        %rotate
        deltaPi = 2*pi/(n/delta);
        curTheta = round((thetaGroup-thetaSubj)/deltaPi);
        for i = 1:prod(size(curHR2,[2 3 4]))
            curHR3(:,i) = circshift(curHR3(:,i),-curTheta);
        end
%         plot(curT3,curHR3(:,i),'-');
        
        %downsample
        curHR4 = nan(size(curHR));
        for i = 1:prod(size(curHR,[2 3 4]))
            curHR4(:,i) = interp1(curT3,curHR3(:,i),t,interpMethod);
        end
%         plot(curT,curHR4(:,i),'o');
        
        hrP{subjInd}.(['sess' num2str(sessInd)]) = permute(curHR4,[2 3 4 1]);
    end
end

% Average runs
condList = {'ori1' 'ori2' 'plaid'};
xData = nan(length(subjList),length(condList),2,12);
for subjInd = 1:length(subjList)
    xData(subjInd,:,:,:) = permute(cat(1,mean(hrP{subjInd}.sess1,1),mean(hrP{subjInd}.sess2,1)),[2 3 1 4]);
end
xDataInfo = 'subj x cond[or1, ori2, plaid] x sess x t';

% Average ori
xData = cat(2,xData,mean(xData(:,1:2,:,:),2));
xDataInfo = 'subj x cond[or1, ori2, plaid, ori] x t';

% Average sessions
xData = squeeze(mean(xData,3));
xDataInfo = 'subj x cond[or1, ori2, plaid, ori] x t';


%% Plot results
fGroup = figure('WindowStyle','docked');
for condInd = [4 3]
    y = squeeze(xData(:,condInd,:));
    t = 0:11;
    yMean = mean(y,1);
    ySEM = std(y,[],1)./sqrt(size(y,1));
    he(condInd) = errorbar(t,yMean,ySEM); hold on
    he(condInd).CapSize = 0;
end
he(4).Color = 'k';
xlabel('time (sec)');
ylabel('?%BOLD change?');
legend(char({'ori +/-SEM' 'plaid +/-SEM'}),'Location','south','box','off','color','w')
title('Group Hemodynamic Response')