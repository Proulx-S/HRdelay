function runGroupAnalysis_hr(saveFig)
close all
if ~exist('saveFig','var') || isempty(saveFig)
    saveFig = 0;
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
% 
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
dCAll = cell(size(subjList,1),1);
for subjInd = 1:size(subjList,1)
    load(fullfile(funPath,funLevel,[subjList{subjInd} fileSuffix]),'dC');
    dCAll{subjInd} = dC;
end
dC = dCAll; clear dAll

% Average voxels
for subjInd = 1:length(subjList)
    for sessInd = 1:2
        % Between-session feature selection
        sess = ['sess' num2str(sessInd)];
        sessFeat = ['sess' num2str(~(sessInd-1)+1)];
        ind = true(1,size(dC{subjInd}.(sess).data,2));
        ind = ind & dC{subjInd}.(sessFeat).anyCondActivation_mask;
        ind = ind & ~dC{subjInd}.(sessFeat).vein_mask;
        
        dC{subjInd}.(sess).data = mean(dC{subjInd}.(sess).data(:,ind,:),2);
        dC{subjInd}.(sess).hr = mean(dC{subjInd}.(sess).hr(:,ind,:,:),2);
    end
end

% Average runs in cartesian space (just to remove between-subject variation in delay from the hr)
condList = {'ori1' 'ori2' 'plaid'};
xData = nan(length(subjList),length(condList),2);
for subjInd = 1:length(subjList)
    xData(subjInd,:,:) = permute(cat(1,mean(dC{subjInd}.sess1.data,1),mean(dC{subjInd}.sess2.data,1)),[2 3 1]);
end
xDataInfo = 'subj x cond[or1, ori2, plaid] x sess';


% Remove subjet effect
rhoGroup = abs(mean(xData(:)));
thetaGroup = angle(mean(xData(:)));
thetaShift = 0; % arbitrary delay (just removing thetaGroup gives phase=0)
xData = mean(xData,2);
t = 0:11;
for subjInd = 1:length(subjList)
    rhoSubj = abs(xData(subjInd,:,:));
    thetaSubj = angle(xData(subjInd,:,:));
    
    interpMethod = 'cubic'; % cubic convolution as implemented in Matlab2020b
    for sessInd = 1:2
        sess = ['sess' num2str(sessInd)];
        % response amplitude
        dC{subjInd}.(sess).hr = dC{subjInd}.(sess).hr./rhoSubj(sessInd).*rhoGroup;
        
        % response delay
        %upsample
        curHR = dC{subjInd}.(sess).hr;
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
        curTheta = round((thetaGroup-thetaSubj(sessInd)+thetaShift)/deltaPi);
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
        
        dC{subjInd}.(sess).hr = permute(curHR4,[2 3 4 1]);
    end
end

% Average runs
condList = {'ori1' 'ori2' 'plaid'};
xData = nan(length(subjList),length(condList),2,12);
for subjInd = 1:length(subjList)
    xData(subjInd,:,:,:) = permute(cat(1,mean(dC{subjInd}.sess1.hr,1),mean(dC{subjInd}.sess2.hr,1)),[2 3 1 4]);
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
xDataPlot = xData*100;
tPlot = t;
xDataPlot = cat(3,xDataPlot,xDataPlot(:,:,1:7));
tPlot = cat(2,tPlot,tPlot(1:7)+n);
% xDataPlot = circshift(xDataPlot,round(thetaShift/(2*pi)*n));
% tPlot = tPlot+round(thetaShift/(2*pi)*n);
clear he
for condInd = [4 3]
    y = squeeze(xDataPlot(:,condInd,:));
    yMean = mean(y,1);
    ySEM = std(y,[],1)./sqrt(size(y,1));
    he(condInd) = errorbar(tPlot,yMean,ySEM); hold on
    he(condInd).CapSize = 0;
end
he(4).Color = 'k';
he(3).Color = colors(3,:);


ht = text(0,0,'ON 6sec');
up = min([he(ishandle(he)).YData]-[he(ishandle(he)).YNegativeDelta]) - diff(ylim)*0.01;
onBarX = [0 0 6 6 0];
down = up-ht.Extent(4)*1.1;
onBarY = [up down down up up];
patch(onBarX,onBarY,'k')
ht.Color = [1 1 1]*0.999999999999;
% ht.Color = [1 1 1];
ht.HorizontalAlignment = 'center'; ht.VerticalAlignment = 'middle';
ht.Position([1 2]) = [mean(onBarX([2 3])) mean(onBarY([1 2]))];
uistack(ht,'top');

onBarX = onBarX+n;
patch(onBarX,onBarY,'k')
ht = copyobj(ht,gca);
ht.Position(1) = mean(onBarX([2 3]));
uistack(ht,'top');

ht = copyobj(ht,gca);
ht.Position(1) = ht.Position(1) - n/2;
ht.Color = 'k';
ht.String = 'OFF 6sec';

BOLDbarX = [0 0 up-down up-down 0];
BOLDbarY = [0.5 0 0 0.5 0.5];
patch(BOLDbarX,BOLDbarY,'k')
text(BOLDbarX(3),mean(BOLDbarY([1 2])),'0.5 %BOLD');

axis off
axis tight
ylabel('?%BOLD change?');
legend(char({'ori +/-SEM' 'plaid +/-SEM'}),'Location','northeast','box','off','color','w')
title('Group Hemodynamic Response')
ax = gca;
uistack(findobj(ax.Children,'type','Text'),'top')
ax.PlotBoxAspectRatio = [1.5 1 1];

if saveFig
    filename = fullfile(pwd,mfilename);
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,'hr');
    fGroup.Color = 'none';
    set(findobj(fGroup.Children,'type','Axes'),'color','none')
    saveas(fGroup,[filename '.svg']); disp([filename '.svg'])
    fGroup.Color = 'w';
    set(findobj(fGroup.Children,'type','Axes'),'color','w')
    saveas(fGroup,[filename '.fig']); disp([filename '.fig'])
    saveas(fGroup,[filename '.jpg']); disp([filename '.jpg'])
end


