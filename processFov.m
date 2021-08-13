function processFov(p,doWhat)
disp('------------')
disp('processFov.m')
p.fov.eccLim = [0.75 7];
outFile = fullfile(p.dataPath.V1,'fov.mat');

switch doWhat
    case {'run' 'run_forced'}
        if ~exist(outFile,'file') || strcmp(doWhat,'run_forced')
            [d,info] = loadData(p);
            disp('empiricalFov: Computing')
            [areaAndFov,cont,voxProp,pEmpirical,fAll,f] = empiricalFov(d,p);
            disp('empiricalFov: Saving')
            save(outFile,'areaAndFov','cont','voxProp','pEmpirical','fAll','f')
            disp(['empiricalFov: Saved to: ' outFile])
        else
            disp('empiricalFov: loading data previously saved locally')
            load(outFile)
        end
    case {'download' 'download_forced'}
        if ~exist(outFile,'file') || strcmp(doWhat,'download_forced')
            disp('empiricalFov: Downloading from repo')
            curLink = 'https://zenodo.org/record/5192849/files/fov.mat?download=1';
            websave(outFile,curLink);
            disp(['empiricalFov: Downloaded to ' outFile])
            disp('empiricalFov: Loading')
        else
            disp('empiricalFov: Loading data previously saved locally')
        end
        load(outFile)
        disp('empiricalFov: Loaded')
end


%% Figure in main paper
figInd = 10;
if p.figOption.verbose>=1
    figInd = [1:11];
end
if p.figOption.save
    figDir = fullfile(p.figOption.outDir,'FigA'); if ~exist(figDir,'dir'); mkdir(figDir); end
    disp('empiricalFov: Saving figures')
end
for i = 1:length(figInd)
    f.L{p.figOption.subjInd,p.figOption.sessInd}(figInd(i)).Visible = 'on';
    f.R{p.figOption.subjInd,p.figOption.sessInd}(figInd(i)).Visible = 'on';
    if p.figOption.save
        figFile = [p.meta.subjList{p.figOption.subjInd} '_L'];
        saveas(f.L{p.figOption.subjInd,p.figOption.sessInd}(figInd(i)),[fullfile(figDir,figFile) '__' num2str(i) '.fig'])
        saveas(f.L{p.figOption.subjInd,p.figOption.sessInd}(figInd(i)),[fullfile(figDir,figFile) '__' num2str(i) '.jpg'])
        figFile = [p.meta.subjList{p.figOption.subjInd} '_R'];
        saveas(f.R{p.figOption.subjInd,p.figOption.sessInd}(figInd(i)),[fullfile(figDir,figFile) '__' num2str(i) '.fig'])
        saveas(f.R{p.figOption.subjInd,p.figOption.sessInd}(figInd(i)),[fullfile(figDir,figFile) '__' num2str(i) '.jpg'])
    end
end
if p.figOption.save
    disp(['empiricalFov: Saved to ' figDir])
end

%% Figure in supplementary
figInd = length(fAll);
if p.figOption.verbose>=1
    figInd = 1:length(fAll);
end
if p.figOption.save
    figDir = fullfile(p.figOption.outDir,'FigSuppl'); if ~exist(figDir,'dir'); mkdir(figDir); end
    disp('empiricalFov: Saving figures')
end
for i = 1:length(figInd)
    fAll{figInd(i)}.Visible = 'on';
    if p.figOption.save
        figFile = [p.meta.subjList{p.figOption.subjInd}];
        saveas(fAll{figInd(i)},[fullfile(figDir,figFile) '__' num2str(i) '.fig'])
        saveas(fAll{figInd(i)},[fullfile(figDir,figFile) '__' num2str(i) '.jpg'])
    end
end
if p.figOption.save
    disp(['empiricalFov: Saved to ' figDir])
end
disp(['empiricalFov: Cleaning up'])
close(findall(groot,'Type','figure','visible','off'))
disp(['empiricalFov: Done'])


function [d,info] = loadData(p)
dataIn = fullfile(p.dataPath.V1,'resp');
dAll = cell(length(p.meta.subjList),1);
for subjInd = 1:length(p.meta.subjList)
    curFile = fullfile(dataIn,[p.meta.subjList{subjInd} '.mat']);
    disp([p.meta.subjList{subjInd} ': loading responses']);
    load(curFile,'resp');
    dAll{subjInd} = resp;
end
disp('responses loaded')
d = dAll; clear dAll
sessList = fields(d{1});
dP = cell(size(d,2),length(sessList));
for subjInd = 1:length(d)
    for sessInd = 1:length(sessList)
        sess = ['sess' num2str(sessInd)];
        dP{subjInd,sessInd} = d{subjInd}.(sessList{sessInd});
        d{subjInd}.(sessList{sessInd}) = [];
    end
end
d = dP; clear dP

info.info = 'subj x sess';
info.sessList = sessList';
info.subjList = p.meta.subjList;
