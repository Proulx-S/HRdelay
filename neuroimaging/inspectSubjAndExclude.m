function inspectSubjAndExclude(figOption,verbose)
close all
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

if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
dataDir = 'C-derived\DecodingHR';
funPath = fullfile(repoPath,dataDir,'fun');
funLevel = 'z';
fileSuffix = '_preprocAndShowMasks.mat';

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



%% Load data
dAll = cell(size(subjList,1),1);
paramAll = cell(size(subjList,1),1);
for subjInd = 1:size(subjList,1)
    curFile = fullfile(funPath,funLevel,[subjList{subjInd} fileSuffix]);
    if verbose; disp(['loading: ' curFile]); end
    load(curFile,'d','param');
    dAll{subjInd} = d;
    paramAll{subjInd} = param;
end


%% Visualize the scanner trigger problem
clear tmp tmp1 tmp2 tmpData
f = figure('WindowStyle','docked');

sz = 0;
for subjInd = 1:length(subjList)
    tmp = size(dAll{subjInd}.sess1.data,1);
    if tmp>sz; sz = tmp; end
end

for subjInd = 1:length(subjList)
    tmp1 = squeeze(mean(dAll{subjInd}.sess1.data,2))';
    tmpInd = squeeze(dAll{subjInd}.sess1.runLabel)';
    [~,b] = sort(tmpInd(:));
    tmpInd(b) = 1:numel(tmpInd);
    tmp1 = tmp1(tmpInd);
    tmp1 = tmp1(:);
    tmp1 = angle(tmp1);
    tmp1 = wrapToPi(tmp1-circ_mean(tmp1));

    tmp2 = squeeze(mean(dAll{subjInd}.sess2.data,2))';
    tmpInd = squeeze(dAll{subjInd}.sess2.runLabel)';
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
dataTmp = angle(mean(dAll{subjInd}.(['sess' num2str(sessInd)]).data,2));
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
        sessInd = exclusion.sess{i}; sess = ['sess' num2str(sessInd)];
        runInd = exclusion.run{i};
        if verbose; disp(['Excluding: ' subjList{subjInd} ', sess' num2str(sessInd) ', runTriplet(repeat)=' num2str(runInd)]); end

        allFields = fields(dAll{subjInd}.(sess));
        nRepeat = size(dAll{subjInd}.(sess).data,1);
        for ii = 1:length(allFields)
            curField = dAll{subjInd}.(sess).(allFields{ii});
            isDataField = ( isnumeric(curField) || islogical(curField) ) && size(curField,1)==nRepeat;
            if isDataField
                dAll{subjInd}.(sess).(allFields{ii})(runInd,:,:,:) = [];
            end
        end
    end
end


%% Plot single subjects
rLim = nan(2,length(subjList),2);
yLim = nan(2,length(subjList),2);
fSubj = cell(length(subjList),1);
for subjInd = 1:length(subjList)
    if subjInd==1
        visibility = 'on';
    elseif figOption.subj==inf && ~verbose
        visibility = 'off';
    else
        visibility = 'on';
    end
    if subjInd==figOption.subj || figOption.subj==inf
        fSubj{subjInd} = figure('WindowStyle','docked','visible',visibility);
        for sessInd = 1:2
            % Between-session feature selection
            sess = ['sess' num2str(sessInd)];
            sessFeat = ['sess' num2str(~(sessInd-1)+1)];
            ind = true(1,size(dAll{subjInd}.(sess).data,2));
            ind = ind & dAll{subjInd}.(sessFeat).anyCondActivation_mask;
            ind = ind & ~dAll{subjInd}.(sessFeat).vein_mask;
            % Polar plot
            subplot(2,2,sessInd)
            xData = squeeze(mean(dAll{subjInd}.(sess).data(:,ind,:),2));
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
            title([subjList{subjInd} '; sess' num2str(sessInd)])
            rLim(:,subjInd,sessInd) = rlim;

            ax = gca;
            ax.ThetaTickLabel = 12-ax.ThetaTick(1:end)/360*12;
            ax.ThetaTickLabel(1,:) = '0 ';

            % Full time course
            subplot(2,2,2+sessInd)
            xData = squeeze(mean(dAll{subjInd}.(sess).hr(:,ind,:,:),2));
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
    end
end
rLim = [0 max(rLim(2,:))];
yLim = [min(yLim(1,:)) max(yLim(2,:))];
for subjInd = 1:length(subjList)
    if subjInd==figOption.subj || figOption.subj==inf
        set(0, 'CurrentFigure', fSubj{subjInd})
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
            fSubj{subjInd}.Color = 'none';
            set(findobj(fSubj{subjInd}.Children,'type','Axes'),'color','none')
            set(findobj(fSubj{subjInd}.Children,'type','PolarAxes'),'color','none')
            curExt = 'svg';
            saveas(fSubj{subjInd},[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
            fSubj{subjInd}.Color = 'w';
            curExt = 'fig';
            saveas(fSubj{subjInd},[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
            curExt = 'jpg';
            saveas(fSubj{subjInd},[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
        end
    else
        if figOption.save
            filename = fullfile(pwd,mfilename);
            if ~exist(filename,'dir'); mkdir(filename); end
            filename = fullfile(filename,subjList{subjInd});
            curFile = dir([filename '.*']);
            if ~isempty(curFile)
                for ii = 1:length(curFile)
                    delFile = fullfile(curFile(ii).folder,curFile(ii).name);
                    if verbose; disp(['delete old: ' delFile]); end
                    delete(delFile)
                end
            end
        end
    end
end


%% Save cleaned data
if verbose; disp('Updating param and cleaned data to:'); end
for subjInd = 1:length(subjList)
    tmp = fullfile(funPath,funLevel,[subjList{subjInd} fileSuffix]);
    if verbose; disp(tmp); end
    dC = dAll{subjInd};
    save(tmp,'dC','param','-append');
end
