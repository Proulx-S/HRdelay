function visualizeFeatSel(p)
pLimAmp = [20 100];
pLimVein = [0 100];
if ~isfield(p,'figOption') || isempty(p.figOption)
    p.figOption.verbose = 1;
    p.figOption.subjInd = 1;
    p.figOption.sessInd = 1;
end
condPair = 'grat1VSgrat2';


%% Define paths
subjList = p.meta.subjList;
repoPath = p.paths.repo.in;
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'd';
            inDir2  = ['e_' p.anaID];
            inDirX = 'c';
            outDir  = ['e_' p.anaID];
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir' 'outDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp


%% Load empiricalFOV figure: for single-subject figure
[fFOV,fDensity] = replotOneFOV(p,funPath,inDir,p.figOption.subjInd,p.figOption.sessInd);

%% Load empiricalFOV figure: for single-subject figure
subjList = p.meta.subjList;
sessList = {'sess1' 'sess2'};
for subjInd = 1:6
    for sessInd = 1:2
        [fFOVall,~] = replotOneFOV(p,funPath,inDir,subjInd,sessInd);
        delete(findobj(fFOVall.Children,'Type','ColorBar'));
        drawnow
        ax = findobj(fFOVall.Children,'Type','Axes');
        for hemiInd = 1:2
            ax(hemiInd).YTick(cellfun('isempty',ax(hemiInd).YTickLabel)) = [];
            ax(hemiInd).YTickLabel(cellfun('isempty',ax(hemiInd).YTickLabel)) = [];
            ax(hemiInd).YTickLabel = [];
            ax(hemiInd).Box = 'off';
            hLine = findobj(ax(hemiInd).Children,'Type','Line');
            set(hLine,'LineWidth',1.5)
            hPoly = findobj(ax(hemiInd).Children,'Type','Polygon');
            set(hPoly,'LineWidth',1.5)
        end
        fFOVall.Color = 'w';
        
        filename = fullfile(p.figOption.finalDir,mfilename);
        filename = fullfile(filename,'FOVs');
        if ~exist(filename,'dir'); mkdir(filename); end
        filename = fullfile(filename,[subjList{subjInd} '_' sessList{sessInd}]);
        
        print(fFOVall,[filename '.png'],'-dpng','-r600')
    end
end


%% Load data
pAll = cell(size(subjList));
for subjInd = 1:length(subjList)
    tmp = load(fullfile(funPath,inDirX,[subjList{subjInd} '.mat']),'p');
    pAll{subjInd} = tmp.p; clear tmp
end
load(fullfile(funPath,inDir2,'featSel.mat'));

%% Plot GLM designs
f = showGLMdesign(featSel{1}.GLMs);
f{1}.Children.DataAspectRatio = [2 10 1];
f{2}.Children.DataAspectRatio = [2 10 1];

f{1}.Children.Units = 'centimeter';
f{2}.Children.Units = 'centimeter';
f{1}.Children.Position(4) = f{2}.Children.Position(4);


% angle(mean(res.(['sess' num2str(p.figOption.sessInd)]).sin(:)))/pi*6
% angle(mean(d.sin(:)))/pi*6

%% Initialize subj, sess and masks
subjInd = p.figOption.subjInd;
sessInd = p.figOption.sessInd;
brain = pAll{subjInd}.brain.mean(:,:,:,sessInd);
roi = pAll{subjInd}.voxProp.area.ind==pAll{subjInd}.voxProp.area.indList(ismember(pAll{subjInd}.voxProp.area.labelList,p.featSel.fov.areaLabel))...
    & ~pAll{subjInd}.voxProp.censorMask;



[ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSel{subjInd,sessInd}.featSeq.featSelList,featSel{subjInd,sessInd}.featSeq.condPairList,p.featSel.global.method,condPair);

featP = nan(size(featSel{subjInd,sessInd}.featSeq.featP,[1 2]));
featP(:,ind_nSpecFeatSel) = featSel{subjInd,sessInd}.featSeq.featP(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond);
featP(:,ind_specFeatSel) = featSel{subjInd,sessInd}.featSeq.featP(:,ind_specFeatSel,ind_specFeatSelCond);

featIndStart = nan(size(featSel{subjInd,sessInd}.featSeq.featIndStart,[1 2]));
featIndStart(:,ind_nSpecFeatSel) = featSel{subjInd,sessInd}.featSeq.featIndStart(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond);
featIndStart(:,ind_specFeatSel) = featSel{subjInd,sessInd}.featSeq.featIndStart(:,ind_specFeatSel,ind_specFeatSelCond);

featQtile = nan(size(featSel{subjInd,sessInd}.featSeq.featQtile,[1 2]));
featQtile(:,ind_nSpecFeatSel) = featSel{subjInd,sessInd}.featSeq.featQtile(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond);
featQtile(:,ind_specFeatSel) = featSel{subjInd,sessInd}.featSeq.featQtile(:,ind_specFeatSel,ind_specFeatSelCond);
featVal = nan(size(featSel{subjInd,sessInd}.featSeq.featVal,[1 2]));
featVal(:,ind_nSpecFeatSel) = featSel{subjInd,sessInd}.featSeq.featVal(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond);
featVal(:,ind_specFeatSel) = featSel{subjInd,sessInd}.featSeq.featVal(:,ind_specFeatSel,ind_specFeatSelCond);
indInSeq = nan(size(featSel{subjInd,sessInd}.featSeq.featIndIn,[1 2]));
indInSeq(:,ind_nSpecFeatSel) = featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond);
indInSeq(:,ind_specFeatSel) = featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_specFeatSel,ind_specFeatSelCond);
indIn = all(indInSeq,2);


featSelIn = false(size(brain));
featSelIn(roi) = indIn;
featSelSeqIn = false([length(featSel{subjInd,sessInd}.featSeq.featSelList) size(brain)]);
featSelSeqIn(:,roi) = indInSeq';
featSelSeqIn = permute(featSelSeqIn,[2 3 4 1]);
featSelVal = nan([size(featSel{subjInd,sessInd}.featSeq.featVal,2) size(brain)]);
featSelValThresholded = nan([size(featSel{subjInd,sessInd}.featSeq.featVal,2) size(brain)]);
featInfo = featSel{subjInd,sessInd}.featSeq.featSelList;
featLabel = cell(size(featInfo));
for featInd = 1:size(featInfo,2)
    tmp = strsplit(featInfo{featInd},': ');
    featLabel(featInd) = tmp(1);
    featSelVal(featInd,roi) = featVal(:,featInd);
    featSelValThresholded(featInd,roi) = featVal(:,featInd);
    featSelValThresholded(featInd,~featSelIn) = nan;
end
featSelVal = permute(featSelVal,[2 3 4 1]);
featSelValThresholded = permute(featSelValThresholded,[2 3 4 1]);





sess = ['sess' num2str(sessInd)];
slice = p.figOption.sliceInd;
visibility = true;

%% Get some nice colormaps
filename = fullfile(pwd,mfilename);
if ~exist(filename,'dir'); mkdir(filename); end
filename = fullfile(filename,'cmap.mat');
if exist(filename,'file')
    load(filename,'cMap_F','cMap_vein');
else
    try
        cMap_F = brewermap(256,'reds');
        cMap_vein = brewermap(256,'blues');
    catch
        error(['Please put this toolbox in .../GitHub/utilities/ (or wherever as long as it ends up in your Matlab path):' newline 'https://github.com/DrosteEffect/BrewerMap'])
    end
    save(filename,'cMap_F','cMap_vein');
end

% %%
% condPairInd = 1;
% % featSel{1}.featSeq.condPairList(condPairInd)
% indIn1 = all(featSel{1}.featSeq.featIndIn(:,1,condPairInd),2);
% indIn2 = all(featSel{1}.featSeq.featIndIn(:,2,condPairInd),2);
% indIn3 = all(featSel{1}.featSeq.featIndIn(:,3,condPairInd),2);
% 
% % nnz(indIn1&indIn2)/nnz(indIn1)
% % nnz(indIn1&indIn3)/nnz(indIn1)
% % 
% % 
% % nnz(indIn1)
% % nnz(indIn2)/nnz(indIn1)
% % nnz(indIn3)
% % nnz(indIn1&indIn2)
% 
% 
% x = featSel{1}.featSeq.featVal(indIn1,2,condPairInd);
% y = featSel{1}.featSeq.featVal(indIn1,3,condPairInd);
% 
% figure('windowstyle','docked')
% hScat = scatter(x,y);
% hScat.MarkerEdgeColor = 'none';
% hScat.MarkerFaceColor = 'k';
% alpha(hScat,0.1)
% ax = gca;
% ax.XAxis.Scale = 'log';
% ax.YAxis.Scale = 'log';


%% Feature selection summary
nVox = nan([5 size(featSel)]);
for i = 1:numel(nVox(1,:))
    nVox(1,i) = size(featSel{i}.featSeq.featIndIn,1);
    nVox(2,i) = nnz(featSel{i}.featSeq.featIndIn(:,1));
    nVox(3,i) = nnz(all(featSel{i}.featSeq.featIndIn(:,1:2),2));
    nVox(4,i) = nnz(all(featSel{i}.featSeq.featIndIn(:,1:3),2));
    nVox(5,i) = nnz(all(featSel{i}.featSeq.featIndIn(:,1:4),2));
end
disp(nVox(:,subjInd,sessInd))
nVox = mean(nVox,3)';
disp(nVox)


%% Plot Brain
f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
axBak = plotIm(axes,brain(:,:,slice));
title('BOLD image (1 TR)')


if 0
    doIt = true(size(featInfo));
    %% Plot unthresholded maps
    for featInd = 1:size(featInfo,2)
        if doIt(featInd)
            f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
            axBak = plotIm(axes,brain(:,:,slice));
            im = featSelVal(:,:,slice,featInd);
            switch featLabel{featInd}
                case 'vein'
                    yLabel = {'Vein score' '(BOLD var / BOLD mean)'};
                    cMap = cMap_vein;
                    scale = 'log';
                case 'act'
                    yLabel = {'Activation' 'Level (F value)'};
                    cMap = cMap_F;
                    scale = 'log';
                case 'respVecSig'
                    yLabel = {'Response Vector Strength' '(Hotelling''s U)'};
                    cMap = cMap_F;
                    scale = 'log';
                case {'respVecDist' 'lateVein'}
                    continue
                case 'respVecDiff'
                    yLabel = {'Response Vector Stimulus Preference' '(Hotelling''s U)'};
                    cMap = cMap_F;
                    scale = 'log';
                case 'retinoFov'
                    yLabel = {'Activation vs Deactivation' '(absolute delay from roi delay)'};
                    cMap = flipud(redblue(256));
                    scale = 'lin';
                    im = im./pi*6;
                otherwise
                    error('X')
            end
            switch scale
                case 'lin'
                case 'log'
                    im = log(im);
                otherwise
                    error('X')
            end
            if strcmp(featLabel{featInd},'retinoFov')
                cLim = [0 6];
            else
                cLim = [min(im(:)) max(im(:))];
                cLim(1) = cLim(1) + diff(cLim)*0.2;
            end
            axOver = plotIm(axes,im,cLim);
            alphaData = ~isnan(im);
            makeOverlay(axBak,axOver,alphaData,cMap,scale,cLim)
            ylabel(yLabel);
            f{end}.Children(1).Position = f{end}.Children(2).Position;
        end
    end
end


%% Plot extra maps
if exist(fullfile(pwd,[mfilename '_crop.mat']),'file')
    load(fullfile(pwd,[mfilename '_crop']),'lim')
else
    lim = [];
end
if ~isfield(lim,['S' p.meta.subjList{p.figOption.subjInd}])
    lim.(['S' p.meta.subjList{p.figOption.subjInd}]) = [];
end
[brainCrop,lim.(['S' p.meta.subjList{p.figOption.subjInd}])] = cropThis(brain,lim.(['S' p.meta.subjList{p.figOption.subjInd}]));
save(fullfile(pwd,[mfilename '_crop']),'lim')

mapInfo = {'none' 'amp' 'vein'};
fExtra = cell(11,1);
indXall = cell(1,size(brain,3));
indYall = cell(1,size(brain,3));
for mapInd = 1:length(mapInfo)    
    for sliceInd = 2:size(brain,3)-1
        if mapInd==1
            fExtra{sliceInd} = figure('WindowStyle','docked','color','w','visible',visibility);
            tExtra{sliceInd} = tiledlayout(length(mapInfo),1,'TileIndexing','columnmajor');
        else
            figure(fExtra{sliceInd})
        end
        curTile = nexttile(tExtra{sliceInd});
        brainSlice = brainCrop(:,:,sliceInd);
        brainSlice(all(brainCrop(:,:,sliceInd)==max(brainCrop(:,:,sliceInd)),2),:) = [];
        brainSlice(:,all(brainCrop(:,:,sliceInd)==max(brainCrop(:,:,sliceInd)),1)) = [];
        
        switch mapInfo{mapInd}
            case 'none'
                axBak = plotIm(curTile,brainSlice);
                continue
            case 'amp'
                tmp = load(fullfile(funPath,inDir,[subjList{p.figOption.subjInd} '.mat']));
                vec = tmp.res.(['sess' num2str(p.figOption.sessInd)]).sin;
                vec = mean(vec(:,:,:,:,p.figOption.condInd),4);
                im = nan(size(roi));
                im(roi) = abs(vec);
                
                tmp = im(:,:,6:end);
                cLim = log(prctile(tmp(:),pLimAmp));
%                 cLim = log([min(im(:)) max(im(:))]);
                im = im(:,:,sliceInd);
                im(all(brainCrop(:,:,sliceInd)==max(brainCrop(:,:,sliceInd)),2),:) = [];
                im(:,all(brainCrop(:,:,sliceInd)==max(brainCrop(:,:,sliceInd)),1)) = [];
                
                im = log(im);
                scale = 'log';
%                 cLim = [min(im(:)) max(im(:))];
%                 if cLim(2)<log(10)
%                     cLim(2) = log(10);
%                 end
                cMap = cMap_F;
                yLabel = {'Response Amplitude' '(%BOLD)'};
            case 'vein'
                im = featSelVal(:,:,:,ismember(featLabel,'vein'));
                tmp = im(:,:,6:end);
                cLim = log(prctile(tmp(:),pLimVein));
%                 cLim = log([min(im(:)) max(im(:))]);
                im = im(:,:,sliceInd);
                im(all(brainCrop(:,:,sliceInd)==max(brainCrop(:,:,sliceInd)),2),:) = [];
                im(:,all(brainCrop(:,:,sliceInd)==max(brainCrop(:,:,sliceInd)),1)) = [];
                im = log(im);
                scale = 'log';
%                 cLim = [min(im(:)) max(im(:))];
                cMap = cMap_vein;
                yLabel = {'Vein score' '(BOLD var / BOLD mean)'};
            otherwise
                error('X')
        end
        
        %Crop out areas that are not in the contiguous ROI
        in = any(~isnan(im),1);
        ind = true(size(in));
        zpos = find(~[0 in 0]);
        [~, grpidx] = max(diff(zpos));
        ind(zpos(grpidx):zpos(grpidx+1)-2) = false;
        indX = ind;
        
        in = any(~isnan(im),2)';
        ind = true(size(in));
        zpos = find(~[0 in 0]);
        [~, grpidx] = max(diff(zpos));
        ind(zpos(grpidx):zpos(grpidx+1)-2) = false;
        indY = ind';
        
%         indX = true(1,size(im,2));
%         indY = true(size(im,1),1);
%         if p.figOption.subjInd==2 && sliceInd==11
%             tmp = find(any(~isnan(im),1),2,'last');
%             indX(find(any(~isnan(im),1),1):tmp(1)) = false;
%             tmp = find(any(~isnan(im),2),2);
%             indY(tmp(2):find(any(~isnan(im),2),1,'last')) = false;
%         else
%             indX(find(any(~isnan(im),1),1):find(any(~isnan(im),1),1,'last')) = false;
%             indY(find(any(~isnan(im),2),1):find(any(~isnan(im),2),1,'last')) = false;
%         end
        brainSlice(indY,:) = [];
        brainSlice(:,indX) = [];
        im(indY,:) = [];
        im(:,indX) = [];
        indYall{sliceInd} = indY;
        indXall{sliceInd} = indX;
        
        axBak = plotIm(curTile,brainSlice);
        axOver = plotIm(axes,im,cLim);
        alphaData = ~isnan(im);
        YTicks = makeOverlay(tExtra{sliceInd}.Children(1),axOver,alphaData,cMap,scale,cLim);
        switch mapInfo{mapInd}
            case {'none' 'amp'}
            case 'vein'
                axOver.YTickLabel = cellstr(num2str(YTicks'));
            otherwise
                error('X')
        end
        ylabel(yLabel);
    end
end
for sliceInd = 2:size(brain,3)-1
    tExtra{sliceInd}.TileSpacing = 'compact';
    for mapInd = 2:length(mapInfo)
        fExtra{sliceInd}.Children(length(mapInfo)-mapInd+1).Position = fExtra{sliceInd}.Children(end).Children(length(mapInfo)-mapInd+1).Position;
    end
    fExtra{sliceInd}.Children(3).Children(3).Children
    tExtra{sliceInd}.Children(3).Visible = 'off';
    fExtra{sliceInd}.Children(3).Visible = 'on';
    
    axes(fExtra{sliceInd}.Children(3).Children(3))
    hold on
    x1 = find(~indXall{sliceInd},1)+0.5;
    x2 = find(~indXall{sliceInd},1,'last')+1.5;
    y1 = find(flip(~indYall{sliceInd}),1)-0.5;
    y2 = find(flip(~indYall{sliceInd}),1,'last')+0.5;
    
    plot([1 1].*x1,[y1 y2],'w')
    plot([1 1].*x2,[y1 y2],'w')
    
    plot([x1 x2],[1 1].*y1,'w')
    plot([x1 x2],[1 1].*y2,'w')
    
    f{end+1} = fExtra{sliceInd};
end



if 0
    %% Plot extra maps 2
    if exist(fullfile(pwd,[mfilename '_crop.mat']),'file')
        load(fullfile(pwd,[mfilename '_crop']),'lim')
    else
        lim = [];
    end
    if ~isfield(lim,['S' p.meta.subjList{p.figOption.subjInd}])
        lim.(['S' p.meta.subjList{p.figOption.subjInd}]) = [];
    end
    [brainCrop,lim.(['S' p.meta.subjList{p.figOption.subjInd}])] = cropThis(brain,lim.(['S' p.meta.subjList{p.figOption.subjInd}]));
    save(fullfile(pwd,[mfilename '_crop']),'lim')
    
    mapInfo = {'none' 'amp' 'vein'};
    fExtra = cell(11,1);
    for mapInd = 1:length(mapInfo)
        for sliceInd = 2:size(brain,3)-1
            if mapInd==1
                fExtra{sliceInd} = figure('WindowStyle','docked','color','w','visible',visibility);
                tExtra{sliceInd} = tiledlayout(length(mapInfo),1,'TileIndexing','columnmajor');
            else
                figure(fExtra{sliceInd})
            end
            curTile = nexttile(tExtra{sliceInd});
            brainSlice = brainCrop(:,:,sliceInd);
            brainSlice(all(brainCrop(:,:,sliceInd)==max(brainCrop(:,:,sliceInd)),2),:) = [];
            brainSlice(:,all(brainCrop(:,:,sliceInd)==max(brainCrop(:,:,sliceInd)),1)) = [];
            
            axBak = plotIm(curTile,brainSlice);
            switch mapInfo{mapInd}
                case 'none'
                    continue
                case 'amp'
                    tmp = load(fullfile(funPath,inDir,[subjList{p.figOption.subjInd} '.mat']));
                    vec = tmp.res.(['sess' num2str(p.figOption.sessInd)]).sin;
                    vec = mean(vec(:,:,:,:,p.figOption.condInd),4);
                    im = nan(size(roi));
                    im(roi) = abs(vec);
                    cLim = [min(im(:)) max(im(:))];
                    im = im(:,:,sliceInd);
                    
                    im(all(brainCrop(:,:,sliceInd)==max(brainCrop(:,:,sliceInd)),2),:) = [];
                    im(:,all(brainCrop(:,:,sliceInd)==max(brainCrop(:,:,sliceInd)),1)) = [];
                    
                    %                 im = log(im);
                    scale = 'lin';
                    %                 cLim = [min(im(:)) max(im(:))];
                    %                 if cLim(2)<log(10)
                    %                     cLim(2) = log(10);
                    %                 end
                    cMap = cMap_F;
                    yLabel = {'Response Amplitude' '(%BOLD)'};
                case 'vein'
                    im = featSelVal(:,:,:,ismember(featLabel,'vein'));
                    cLim = [min(im(:)) max(im(:))];
                    im = im(:,:,sliceInd);
                    im(all(brainCrop(:,:,sliceInd)==max(brainCrop(:,:,sliceInd)),2),:) = [];
                    im(:,all(brainCrop(:,:,sliceInd)==max(brainCrop(:,:,sliceInd)),1)) = [];
                    %                 im = log(im);
                    scale = 'lin';
                    %                 cLim = [min(im(:)) max(im(:))];
                    cMap = cMap_vein;
                    yLabel = {'Vein score' '(BOLD var / BOLD mean)'};
                otherwise
                    error('X')
            end
            axOver = plotIm(axes,im,cLim);
            alphaData = ~isnan(im);
            
            makeOverlay(tExtra{sliceInd}.Children(1),axOver,alphaData,cMap,scale,cLim)
            ylabel(yLabel);
        end
    end
    for sliceInd = 2:size(brain,3)-1
        tExtra{sliceInd}.TileSpacing = 'compact';
        for mapInd = 2:length(mapInfo)
            fExtra{sliceInd}.Children(length(mapInfo)-mapInd+1).Position = fExtra{sliceInd}.Children(end).Children(length(mapInfo)-mapInd+1).Position;
        end
            f{end+1} = fExtra{sliceInd};
    end
end



if 0
    %% Plot thresholded maps
    for featInd = 1:size(featInfo,2)
        if doIt(featInd)
            f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
            axBak = plotIm(axes,brain(:,:,slice));
            im = featSelValThresholded(:,:,slice,featInd);
            switch featLabel{featInd}
                case 'vein'
                    yLabel = {'Thresholded Vein score' 'BOLD var / BOLD mean'};
                    cMap = cMap_vein;
                    scale = 'log';
                case 'act'
                    yLabel = {'Thresholded Activation Level' '(F value)'};
                    cMap = cMap_F;
                    scale = 'log';
                case 'respVecSig'
                    yLabel = {'Thresholded Response Vector Strength' '(Hotelling''s U)'};
                    cMap = cMap_F;
                    scale = 'log';
                case {'respVecDist' 'lateVein'}
                    continue
                case 'respVecDiff'
                    yLabel = {'Thresholded Discriminativeness' '(Hotelling''s U)'};
                    cMap = cMap_F;
                case 'retinoFov'
                    yLabel = {'Thresholded Activation vs Deactivation' '(absolute delay from roi delay)'};
                    cMap = flipud(redblue(256));
                    scale = 'lin';
                    im = im./pi*6;
                otherwise
                    error('X')
            end
            switch scale
                case 'lin'
                case 'log'
                    im = log(im);
                otherwise
                    error('X')
            end
            if strcmp(featLabel{featInd},'retinoFov')
                cLim = [0 6];
            else
                cLim = [min(im(:)) max(im(:))];
                cLim(1) = cLim(1) + diff(cLim)*0.2;
            end
            axOver = plotIm(axes,im,cLim);
            alphaData = ~isnan(im);
            makeOverlay(axBak,axOver,alphaData,cMap,scale,cLim)
            ylabel(yLabel)
        end
    end
end

if 0
    %% Plot response vector distribution
    tmp = load(fullfile(funPath,'d',[subjList{subjInd} '.mat']));
    d = tmp.res.(sess); clear tmp
    f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
    x = mean(d.sin(:,:,:,:),4);
    featInd = ismember(featLabel,'respVecSig');
    ind = featSel{subjInd,sessInd}.featSeq.featIndIn(:,featInd,condPairInd);
    polarscatter(angle(x(ind)),abs(x(ind)),eps,'.k'); hold on
    polarscatter(angle(x(~ind)),abs(x(~ind)),eps,'.r');
    
    %% Plot response vector maps
    p.dataType = 'sin';
    p.svmSpace = 'cartRoi';
    voxFlag = 1;
    f0 = plotNorm(d,p,featSel{subjInd,sessInd},voxFlag);
    
    % for slice = 2:12
    f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
    axBak = plotIm(axes,brain(:,:,slice));
    [x,y,~] = getXYK(d,p);
    [x,~] = polarSpaceNormalization(x,p.svmSpace);
    x = real(x);
    x = mean(x(y==2,:),1) - mean(x(y==1,:),1);
    im = nan(size(brain));
    im(roi) = x;
    im = im(:,:,slice);
    yLabel = {'%BOLD amplitude at roi delay'};
    cMap = jet;
    scale = 'lin';
    switch scale
        case 'lin'
        case 'log'
            im = log(im);
        otherwise
            error('X')
    end
    cLim = [-1 1].*max(abs(x(:)));
    cLim = [cLim(1) + diff(cLim)*0.2 cLim(2) - diff(cLim)*0.2];
    axOver = plotIm(axes,im,cLim);
    alphaData = ~isnan(im);
    makeOverlay(axBak,axOver,alphaData,cMap,scale,cLim)
    ylabel([yLabel {[num2str(cLim(1)) ' to ' num2str(cLim(2))]}])
    % end
end

%% Plot histograms
if 0
    doIt = true(size(featInfo));
    fovSelIn = featSelSeqIn(:,:,:,ismember(featLabel,'retinoFov'));
    for featInd = 1:size(featInfo,2)
        if doIt(featInd)
            if strcmp(featLabel{featInd},'retinoFov')
                thresh = p.featSel.fov;
            else
                thresh = p.featSel.(featLabel{featInd});
            end
            switch featLabel{featInd}
                case 'retinoFov'
                    xLabel = {'Delay relative to roi average' '(sec)'};
                    qtile = nan;
                    scale = 'lin';
                    xLim = 'auto';
                case 'vein'
                    xLabel = {'Vein score' 'log( BOLD var / BOLD mean )'};
                    qtile = 1-thresh.percentile/100;
                    scale = 'log';
                    xLim = 'auto';
                case 'act'
                    xLabel = {'Activation Level' 'log( F value )'};
                    qtile = thresh.percentile/100;
                    scale = 'log';
                    xLim = 'auto';
                case 'respVecSig'
                    xLabel = {'Response Vector Strenght' 'log( Hotelling''s U )'};
                    qtile = thresh.percentile/100;
                    scale = 'log';
                    xLim = 'auto';
                case 'lateVein'
                    xLabel = {'Response Vector Delay' '(sec)'};
                    qtile = 1 - thresh.percentile/100;
                    scale = 'lin';
                    xLim = 'auto';
                case 'respVecDist'
                    continue
                case 'respVecDiff'
                    xLabel = {'Discriminativeness' 'log( Hotelling''s U )'};
                    qtile = thresh.percentile/100;
                    scale = 'log';
                    xLim = 'auto';
                otherwise
                    error('X')
            end
            f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
            tmp = featSelVal(:,:,:,featInd);
            if strcmp(featLabel{featInd},'retinoFov')
                tmp = tmp./pi*6;
            end
            switch scale
                case 'log'
                    tmp = log(tmp);
                case 'lin'
                otherwise
                    error('X')
            end
            
            [~,edges] = histcounts(tmp);
            Nin = histcounts(tmp(roi & featSelIn),edges)';
            Nout_fovIn = histcounts(tmp(roi & ~featSelIn & fovSelIn),edges)';
            Nout_fovOut = histcounts(tmp(roi & ~featSelIn & ~fovSelIn),edges)';
            ctrs = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
            hBar = bar(ctrs, [Nin Nout_fovIn Nout_fovOut],1,'stacked','FaceAlpha',0.6); hold on
            hBar(1).FaceColor = 'w';
            hBar(1).EdgeColor = 'none';
            hBar(2).FaceColor = 'r';
            hBar(2).EdgeColor = 'none';
            hBar(3).FaceColor = 'b';
            hBar(3).EdgeColor = 'none';
            
            outline_x = ctrs-mode(diff(edges)/2);
            outline_y = sum([Nin Nout_fovIn Nout_fovOut],2)';
            if outline_y(1)~=0
                outline_x = [ctrs(1)-mode(diff(edges)/2)-mode(diff(edges)) outline_x];
                outline_y = [0 outline_y];
            end
            if outline_y(end)~=0
                outline_x = [outline_x ctrs(end)-mode(diff(edges)/2)+mode(diff(edges))];
                outline_y = [outline_y 0];
            end
            stairs(outline_x,outline_y,'k')
            
            ax = gca;
            ax.Box = 'off';
            
            if ~strcmp(thresh.threshMethod,'empirical')
                featVals = featVal(:,featInd);
                qtiles = featQtile(:,featInd);
                ps = featP(:,featInd);
                fdrs = nan(size(ps));
                startInds = logical(featIndStart(:,featInd));
                switch thresh.threshMethod
                    case '%ile'
                        featValThresh = interp1(qtiles(~isnan(qtiles)&~isnan(featVals)),featVals(~isnan(qtiles)&~isnan(featVals)),qtile);
                    case 'p'
                        featValThresh = interp1(ps,featVals,thresh.threshVal);
                    case 'fdr'
                        fdrs(startInds) = mafdr(ps(startInds),'BHFDR',true);
                        [~,b] = sort(abs(fdrs-thresh.threshVal));
                        featValThresh = mean(featVals(b(1:2)));
                    otherwise
                        error('X')
                end
                switch scale
                    case 'log'
                        featValThresh = log(featValThresh);
                    case 'lin'
                end
                hThresh = plot([1 1].*featValThresh,ylim,':k');
                switch thresh.threshMethod
                    case {'p' 'fdr'}
                        legend([hBar hThresh],{'incl.' 'excl. (within fov)' 'excl. (outside fov)' [thresh.threshMethod '=' num2str(thresh.threshVal)]},'box','off');
                    case '%ile'
                        legend([hBar hThresh],{'incl.' 'excl. (within fov)' 'excl. (outside fov)' [thresh.threshMethod '=' num2str(thresh.percentile)]},'box','off');
                    otherwise
                        error('X')
                end
            else
                legend(hBar,{'incl.' 'excl. (within fov)' 'excl. (outside fov)'},'box','off');
            end
            ax.TickDir = 'out';
            xlim(xLim)
            xlabel(xLabel)
        end
    end
end

if 0
    figure('WindowStyle','docked');
    corrplot(featVal(indIn,:),'varnames',featLabel)
end

%% Save figures
f{end+1} = fFOV;
f{end+1} = fDensity;
f{end+1} = fFOVall;
if 1
    filename = fullfile(p.figOption.finalDir,mfilename);
    
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,[subjList{subjInd}]);
    for i = 1:length(f)
        if i == length(f)
            keyboard
        end
        curF = f{i};
%         curF.Color = 'none';
%         set(findobj(curF.Children,'type','Axes'),'color','none')
        curFile = [filename '_' num2str(i)];
        curExt = 'svg';
        saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
        curExt = 'eps';
        exportgraphics(curF,[curFile '.' curExt],'ContentType','vector')
%         curF.Color = 'w';
        curExt = 'fig';
        saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
        curExt = 'jpg';
        saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
    end
end
