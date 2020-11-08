function curData = getWholeROIdata(subj,subjStim,p)
%% Set some paths and names
mainDir = p.mainDir;
labelDir = p.labelDir;

maskName = 'v1';%v1, v2 or v4
dataType = 'trial'; % 'trial' 'run'
fitType = 'crossCorr';% 'crossCorr' 'corr' 'corrFixDelay' 'corrVoxSpecificDelay'
dataStringIn1 = 'run2_crossCor__';
dataStringIn2 = '';
dataDirIn = fullfile(mainDir,subj);
dataFileIn = [dataStringIn1 subj '_' maskName dataStringIn2];
dataStringOut = 'run4_wholeROIdist';

plotIt = 1;
% doPerm = 0;


p.thresholdData = 1;
p.PL = 1;% phase-locking thresholding

p.trialsPerObs = 8;
p.thresholdP = 0.05;
p.maskName = maskName;
p.dataType = [dataType num2str(p.trialsPerObs)];
p.fitType = fitType;
p.dataFileIn = fullfile(dataDirIn,dataFileIn);
%Should set dataFileOut better
titleStr = p.hist.data;
titleStr = [titleStr '_' 'avGrat' num2str(p.hist.averageGrat)];
% titleStr = [titleStr '_' p.hist.distType];
titleStr = [titleStr '_' 'avRun' num2str(p.hist.averageRun)];
titleStr = [titleStr '__' subj];
 
p.dataFileOut = fullfile(mainDir,[dataStringOut '__' titleStr]);

%% Get data
load(p.dataFileIn,'metric')

tmp.run = metric.run;
tmp.runList = metric.runList;
if ~strcmp(p.dataType,'run')
    tmp.(p.dataType) = metric.(p.dataType);
end
metric = tmp; clear tmp


%Compile it
runList = dir(fullfile(labelDir,subjStim,'/*.mat'));
label = [];
i=1;
while i <= length(runList)
    tmp = strsplit(runList(i).name,'__');
    if length(tmp)==3
        label = [label str2num(tmp{3}(18:end-4))];
    end
    i = i+1;
end

if strcmp(p.dataType,'run')
    metricCompiled = compileRuns(metric,label,p.fitType,p.dataType);
    metricCompiled_thresh = metricCompiled;
else
    metricCompiled = compileRuns(metric,label,p.fitType,p.dataType);
    metricCompiled_thresh = compileRuns(metric,label,p.fitType,'run');    
end

%Threshold it
%based on r^2
if p.thresholdData
    switch p.thresholdP
        case 0.05
            threshRsquared = 0.0358;
        case 0.01
            threshRsquared = 0.0610;
        case 0.001
            threshRsquared = 0.0975;
        case 0.0001
            threshRsquared = 0.1336;
        otherwise
            warning('threshRsquared not precompiled for requested p, will estimate it. Should be fine as long as you have very large noisy data set')
            threshRsquared = findRsqaredAtP(metricCompiled,p.thresholdP);
    end
    metricCompiled = applyThresh(metricCompiled,metricCompiled_thresh,'Rsquared',threshRsquared);%0.2435^2
    metric = metricCompiled; clear metricCompiled
else
    clear metricCompiled
end
%based on phase locking
if p.PL
    metric = applyPhaseLockingThresh(metric);
end

%% Show whole-ROI histograms
p.hist.nBins = 30;
switch p.hist.data
    case 'delayRect'
        p.hist.periodShift = 0.25;
    otherwise
        p.hist.periodShift = 0;
end

curData = metric.(p.hist.data);
%Temporal shift
if p.hist.periodShift
    curData = curData+metric.fittingParam.period*0.25;
end
%Handle runs
if p.hist.averageRun
    switch p.hist.data
        case 'ampRect'
            curData = mean(curData,2);
        case 'delayRect'
            curData = sec2rad(curData,metric.fittingParam.period);
            curData = circ_mean(curData,[],2);
            curData = rad2sec(curData,metric.fittingParam.period);
    end
end
curData = squeeze(curData);

