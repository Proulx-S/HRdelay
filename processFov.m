function processFov(p,doWhat)
disp('------------')
disp('processFov.m')
p.fov.eccLim = [0.75 7];
outFile = fullfile(p.dataPath.V1,'fov.mat');





% [d,info] = loadData(p);
% load(outFile);

% fHist = figure;
% y = angle(mean(d{p.figOption.subjInd,p.figOption.sessInd}.sin(:,:),2));
% histogram(y);
% hold on

% y = [cont.L{p.figOption.subjInd,p.figOption.sessInd}.vecUV; cont.R{p.figOption.subjInd,p.figOption.sessInd}.vecUV];
% histogram(y);


% y = angle(mean(d{p.figOption.subjInd,p.figOption.sessInd}.sin(:,:),2));
% histogram(y);

% fHist = figure;
% histogram(y);
% hold on;
% y = areaAndFov{p.figOption.subjInd,p.figOption.sessInd}.featVal;
% i = logical(areaAndFov{p.figOption.subjInd,p.figOption.sessInd}.featIndIn);
% nnz(i)
% length(y)
% length(cont.L{p.figOption.subjInd,p.figOption.sessInd}.vecUV)+length(cont.R{p.figOption.subjInd,p.figOption.sessInd}.vecUV)
% histogram(y(i));
% histogram(cont.L{p.figOption.subjInd,p.figOption.sessInd}.vecUV);






[d,info] = loadData(p);
switch doWhat
    case {'run' 'run_forced' 'run_forced_butSkipFlattenEccDist'}
        if ~exist(outFile,'file') || strcmp(doWhat,'run_forced') || strcmp(doWhat,'run_forced_butSkipFlattenEccDist')
            disp('empiricalFov: Computing')
            if strcmp(doWhat,'run_forced_butSkipFlattenEccDist')
                load(outFile,'voxProp');
                [areaAndFov,cont,voxProp,pEmpirical,fAll,f] = empiricalFov(d,p,voxProp);
            else
                [areaAndFov,cont,voxProp,pEmpirical,fAll,f] = empiricalFov(d,p);
            end
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


%% Figures
% in main paper
figDir = fullfile(p.figOption.outDir); if ~exist(figDir,'dir'); mkdir(figDir); end
disp('empiricalFov: Saving figures')
figInd = 10;
f.L{p.figOption.subjInd,p.figOption.sessInd}(figInd).Visible = 'on';
f.R{p.figOption.subjInd,p.figOption.sessInd}(figInd).Visible = 'on';
figFile = 'Fig3Aleft';
saveas(f.L{p.figOption.subjInd,p.figOption.sessInd}(figInd),[fullfile(figDir,figFile) '.fig'])
saveas(f.L{p.figOption.subjInd,p.figOption.sessInd}(figInd),[fullfile(figDir,figFile) '.jpg'])
figFile = 'Fig3Aright';
saveas(f.R{p.figOption.subjInd,p.figOption.sessInd}(figInd),[fullfile(figDir,figFile) '.fig'])
saveas(f.R{p.figOption.subjInd,p.figOption.sessInd}(figInd),[fullfile(figDir,figFile) '.jpg'])
% phase histogram inset for Fig3A
phaseDelay = areaAndFov{p.figOption.subjInd,p.figOption.sessInd}.phaseDelay - areaAndFov{p.figOption.subjInd,p.figOption.sessInd}.phaseDelayNorm;
phaseDelay = wrapToPi(-phaseDelay-pi/2)+pi/2;
f.hist = figure('Visible','off');
histogram(phaseDelay,linspace(-pi/2,pi+pi/2,50)); axis tight
xline([0 pi]);
xlabel('phase delay relative to phase of average response (rad)');
figFile = 'Fig3Ahist';
f.hist.Visible = 'on';
saveas(f.hist,[fullfile(figDir,figFile) '.fig'])
saveas(f.hist,[fullfile(figDir,figFile) '.jpg'])



% in supplement
fAll{end}.Visible = 'on';
figFile = 'SuppFig1';
saveas(fAll{end},[fullfile(figDir,figFile) '.fig'])
saveas(fAll{end},[fullfile(figDir,figFile) '.jpg'])
% phase histogram inset for SuppFig1
fAllHist = figure('Visible','off'); ax = cell(size(areaAndFov,1),size(areaAndFov,2));
clear dd ddPhase
for subjInd = 1:size(areaAndFov,1)
    for sessInd = 1:size(areaAndFov,2)
        dd{subjInd,sessInd} = mean(d{subjInd,sessInd}.sin(:,:),2);
        dd{subjInd,sessInd} = dd{subjInd,sessInd}./abs(mean(dd{subjInd,sessInd}));
        ddPhase{subjInd,sessInd} = angle(mean(dd{subjInd,sessInd}));
        dd{subjInd,sessInd} = dd{subjInd,sessInd}.*exp(-1i*ddPhase{subjInd,sessInd});
        ax{subjInd,sessInd} = subplot(size(areaAndFov,2),size(areaAndFov,1),subjInd+(sessInd-1)*size(areaAndFov,1));
        % % phaseDelay{subjInd,sessInd} = areaAndFov{subjInd,sessInd}.phaseDelay - areaAndFov{subjInd,sessInd}.phaseDelayNorm;
        % % phaseDelay{subjInd,sessInd} = wrapToPi(phaseDelay{subjInd,sessInd}-pi/2)+pi/2;
        histogram(wrapToPi(-angle(dd{subjInd,sessInd})-pi/2)+pi/2,linspace(-pi/2,pi+pi/2,50));
        xline([0 pi]);
        axis tight
    end
    drawnow
    yLim = get([ax{subjInd,:}],'YLim');
    yLim = [0 max([yLim{:}])]
    set([ax{subjInd,:}],'YLim',yLim)
    xlabel('vox phase - roi phase (rad)')
end
figFile = 'SuppFig1hist';
fAll{end}.Visible = 'on';
saveas(fAllHist,[fullfile(figDir,figFile) '.fig'])
saveas(fAllHist,[fullfile(figDir,figFile) '.jpg'])

% group phase histogram for fun
for subjInd = 1:size(areaAndFov,1)
    dd{subjInd,1}      = mean(cat(2,dd{subjInd,:}     ),2);
    ddPhase{subjInd,1} = mean(cat(2,ddPhase{subjInd,:}),2);
    % dd{subjInd,1}      = dd{subjInd,1}.*exp(1i*ddPhase{subjInd,1});
end
dd = cat(1,dd{:,1});
% polarhistogram(angle(dd),100);
fCatHist = figure('Visible','off');
histogram(wrapToPi(angle(dd)-pi/2)+pi/2,linspace(-pi/2,pi+pi/2,50)); axis tight
xline([0 pi]);
figFile = 'groupPhaseHist';
saveas(fCatHist,[fullfile(figDir,figFile) '.fig'])
saveas(fCatHist,[fullfile(figDir,figFile) '.jpg'])

disp(['empiricalFov: Saved to ' figDir])
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
