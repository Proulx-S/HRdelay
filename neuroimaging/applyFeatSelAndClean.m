function applyFeatSelAndClean(figOption)
close all
if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 1;
    figOption.subj = 1; % 'all' or subjInd
end

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
fileSuffix = '_defineAndShowMasks.mat';

%make sure everything is forward slash for mac, linux pc compatibility
for tmpPath = {'repoPath' 'dataDir' 'funPath'}
    eval([char(tmpPath) '(strfind(' char(tmpPath) ',''\''))=''/'';']);
end

%% Preload param
tmp = dir(fullfile(funPath,funLevel,'*.mat'));
for i = 1:length(tmp)
    load(fullfile(funPath,funLevel,tmp(i).name),'param')
    if exist('param','var')
        break
    end
end
if ~exist('param','var')
    error('Analysis parameters not found!')
end
subjList = param.subjList;
featSelContrast1 = param.featSelContrast1.name;
threshType = param.featSelContrast1.threshType;
threshVal = param.featSelContrast1.threshVal;


disp(['IN: Sinusoidal BOLD responses from anatomical V1 ROI (' fullfile(dataDir,funLevel) ')'])
disp(['threshVal=' num2str(threshVal)])
disp('F(IN)=OUT: threshold included voxels and analyse responses averaged across the ROI')
disp(['OUT: figures and stats'])



%% Load data
dAll = cell(size(subjList,1),1);
featSelStatsAll = cell(size(subjList,1),1);
veinAll = cell(size(subjList,1),1);
for subjInd = 1:size(subjList,1)
    load(fullfile(funPath,funLevel,[subjList{subjInd} fileSuffix]),'d','featSelStats','vein');
    dAll{subjInd} = d;
    featSelStatsAll{subjInd} = featSelStats;
    veinAll{subjInd} = vein;    
end
d = dAll; clear dAll
featSelStats = featSelStatsAll; clear featSelStatsAll
vein = veinAll; clear veinAll


save('C:\Users\sebas\Desktop\main','d','featSelStats','vein')
return


dC = d;
%% Apply feature selection (cross-session)
rLim = nan(length(dC),1);
for subjInd = 1:length(dC)
    for sessInd = 1:2
        sess = ['sess' num2str(sessInd)];
        sessFeat = ['sess' num2str(~(sessInd-1)+1)];
        % remove veins
        indVein = vein{subjInd}.(sessFeat).mask;
        % remove inactive voxels
        if ~strcmp(threshType,'none')
            indAct = featSelStats{subjInd}.(featSelContrast1).(sessFeat).(upper(threshType))<threshVal;
        else
            indAct = true(size(featSelStats{subjInd}.(featSelContrast1).(sessFeat).(upper(threshType))));
        end
        dC{subjInd}.(sess).data = dC{subjInd}.(sess).data(:,indAct&~indVein,:,:);
        dC{subjInd}.(sess).hr = dC{subjInd}.(sess).hr(:,indAct&~indVein,:,:);
        switch featSelContrast1
            case 'anyCondActivation'
                dC{subjInd}.(sess).([featSelContrast1 '__F']) = featSelStats{subjInd}.(featSelContrast1).(sessFeat).F(:,indAct&~indVein);
            otherwise
                error('X')
        end
    end
    disp('+++')
    if all(~indVein)
        disp('NOT masking out veins vox')
    else
        disp('masking out veins vox')
    end
    if all(indAct)
        disp('NOT masking out inactive vox')
    else
        disp('masking out inactive vox')
    end
    disp(['Leaving ' num2str(size(dC{subjInd}.sess1.data,2)) ' and ' num2str(size(dC{subjInd}.sess2.data,2)) ' vox in sess1 and sess2']);
    disp('+++')
end

%% Exclude the scanner trigger problem
clear tmp tmp1 tmp2 tmpData
f = figure('WindowStyle','docked','visible','off');

sz = 0;
for subjInd = 1:length(subjList)
    tmp = size(dC{subjInd}.sess1.data,1);
    if tmp>sz; sz = tmp; end
end

for subjInd = 1:length(subjList)
    tmp1 = squeeze(mean(dC{subjInd}.sess1.data,2))';
    tmpInd = squeeze(dC{subjInd}.sess1.runLabel)';
    [~,b] = sort(tmpInd(:));
    tmpInd(b) = 1:numel(tmpInd);
    tmp1 = tmp1(tmpInd);
    tmp1 = tmp1(:);
    tmp1 = angle(tmp1);
    tmp1 = wrapToPi(tmp1-circ_mean(tmp1));
    
    tmp2 = squeeze(mean(dC{subjInd}.sess2.data,2))';
    tmpInd = squeeze(dC{subjInd}.sess2.runLabel)';
    [~,b] = sort(tmpInd(:));
    tmpInd(b) = 1:numel(tmpInd);
    tmp2 = tmp2(tmpInd);
    tmp2 = tmp2(:);
    tmp2 = angle(tmp2);
    tmp2 = wrapToPi(tmp2-circ_mean(tmp2));
    
    tmp = cat(1,tmp1,nan(sz*3-length(tmp1)+1,1),tmp2);
    
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

% Update exclude
subjInd = 2;
sessInd = 1;
dataTmp = angle(mean(dC{subjInd}.(['sess' num2str(sessInd)]).data,2));
dataTmp = wrapToPi(dataTmp - circ_mean(dataTmp));
dataTmp = abs(dataTmp);
a = max(dataTmp(:),[],1);
[runInd,condInd] = find(squeeze(dataTmp)==a);
exclusion.subjList = subjList;
exclusion.subj = subjInd;
exclusion.sess = {sessInd};
exclusion.run = {runInd};
exclusion.cond = {1:3};
param.exclusion = exclusion;

% exclusion.subjList
% exclusion.subj = 2;
% exclusion.sess = {1};
% exclusion.run = {4};
% exclusion.cond = {1 2 3};


%% Apply exclusion
if exist('exclusion','var') && ~isempty(exclusion) && ~isempty(exclusion.subj)
    if length(exclusion.subj)>1; error('Not coded for multiple exclusions'); end
    for i = 1:length(exclusion.subj)
        subjInd = find(ismember(subjList,exclusion.subjList{exclusion.subj(i)}));
        sessInd = exclusion.sess{i};
        runInd = exclusion.run{i};
        disp(['Excluding: ' subjList{subjInd} ', sess' num2str(sessInd) ', runTriplet(repeat)=' num2str(runInd)])
        
        dC{subjInd}.(['sess' num2str(sessInd)]).data(runInd,:,:) = [];
        dC{subjInd}.(['sess' num2str(sessInd)]).hr(runInd,:,:,:) = [];
        dC{subjInd}.(['sess' num2str(sessInd)]).runLabel(runInd,:,:) = [];
    end
end


%% Plot single subjects
rLim = nan(2,length(subjList),2);
yLim = nan(2,length(subjList),2);
for subjInd = 1:length(subjList)
    if subjInd==figOption.subj || figOption.subj==inf
        fSubj(subjInd) = figure('WindowStyle','docked');
        for sessInd = 1:2
            % Polar plot
            subplot(2,2,sessInd)
            xData = squeeze(mean(dC{subjInd}.(['sess' num2str(sessInd)]).data,2));
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
            rLim(:,subjInd,sessInd) = rlim;
            
            ax = gca;
            ax.ThetaTickLabel = 12-ax.ThetaTick(1:end)/360*12;
            ax.ThetaTickLabel(1,:) = '0 ';
            
            % Full time course
            subplot(2,2,2+sessInd)
            xData = squeeze(mean(dC{subjInd}.(['sess' num2str(sessInd)]).hr,2));
            hPlot = [];
            hPlotMean = [];
            for condInd = 1:size(xData,2)
                y = squeeze(xData(:,condInd,:));
                t = 0:11;
%                 yMean = mean(y,1);
%                 ySEM = std(y,[],1)./sqrt(size(y,1));
%                 he = errorbar(t,yMean,ySEM); hold on
%                 he.CapSize = 0;
%                 he.Color = colors(condInd,:);
                hPlot(:,condInd) = plot(t,y,':','color',colors(condInd,:)); hold on
                hPlotMean(:,condInd) = plot(t,mean(y,1),'color',colors(condInd,:));
            end
            uistack(hPlotMean,'top')
            set(hPlotMean,'linewidth',1.5)
            box off
            yLim(:,sessInd,subjInd) = ylim;
        end
        
        suptitle(subjList{subjInd})
    else
    end
end
rLim = [0 max(rLim(2,:))];
yLim = [min(yLim(1,:)) max(yLim(2,:))];
for subjInd = 1:length(subjList)
    if subjInd==figOption.subj || figOption.subj==inf
        figure(fSubj(subjInd));
        for sessInd = 1:2
            ax = subplot(2,2,sessInd);
            rlim(rLim);
            switch sessInd
                case 1
                    ax.ThetaAxis.Label.String = {'delay' '(sec)'};
                    ax.ThetaAxis.Label.FontSize = 9;
                    ax.ThetaAxis.Label.Position(1) = sum(ax.Position([1 3]));
                    ax.RAxis.Label.String = '%BOLD';
                    ax.RAxis.Label.Rotation = 80;
                    ax.RAxis.Label.FontSize = 9;
                case 2
                    hl = legend(char({'ori1' 'ori2' 'plaid'}),'Location','east','Color','none','Box','off');
            end
            ax = subplot(2,2,2+sessInd);
            switch sessInd
                case 1
                    ylabel('%BOLD')
                    xlabel('time since stim onset (sec)')
                case 2
                    xlabel('time since stim onset (sec)')
            end
            
            ylim(yLim);
        end
        if figOption.save
            filename = fullfile(pwd,mfilename);
            if ~exist(filename,'dir'); mkdir(filename); end
            filename = fullfile(filename,subjList{subjInd});
            curFile = filename;
            fSubj(subjInd).Color = 'none';
            set(findobj(fSubj(subjInd).Children,'type','Axes'),'color','none')
            set(findobj(fSubj(subjInd).Children,'type','PolarAxes'),'color','none')
            curExt = 'svg';
            saveas(fSubj(subjInd),[curFile '.' curExt]); disp([curFile '.' curExt])
            fSubj(subjInd).Color = 'w';
            curExt = 'fig';
            saveas(fSubj(subjInd),[curFile '.' curExt]); disp([curFile '.' curExt])
            curExt = 'jpg';
            saveas(fSubj(subjInd),[curFile '.' curExt]); disp([curFile '.' curExt])
        end
    else
        if figOption.save
            filename = fullfile(pwd,mfilename);
            if ~exist(filename,'dir'); mkdir(filename); end
            filename = fullfile(filename,subjList{subjInd});
            curFile = filename;
            curExt = 'svg';
            delete([curFile '.' curExt]); disp(['del old: ' curFile '.' curExt])
            curExt = 'fig';
            delete([curFile '.' curExt]); disp(['del old: ' curFile '.' curExt])
            curExt = 'jpg';
            delete([curFile '.' curExt]); disp(['del old: ' curFile '.' curExt])
        end
    end
end


%% Save cleaned data
dCAll = dC;
disp('Updating param and cleaned data to:')
for subjInd = 1:length(subjList)
    tmp = fullfile(funPath,funLevel,[subjList{subjInd} fileSuffix]);
    disp(tmp)
    dC = dCAll{subjInd};
    save(tmp,'dC','param','-append');
end

