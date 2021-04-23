function visualizeFeatSel(p)
if ~isfield(p,'figOption') || isempty(p.figOption)
    p.figOption.verbose = 1;
    p.figOption.subjInd = 1;
    p.figOption.sessInd = 1;
end
condPair = 'grat1VSgrat2';


%% Define paths
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'd';
            inDirX = 'c';
            outDir  = 'd';
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir' 'outDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp

%% Load data
pAll = cell(size(subjList));
for subjInd = 1:length(subjList)
    tmp = load(fullfile(funPath,inDirX,[subjList{subjInd} '.mat']),'p');
    pAll{subjInd} = tmp.p; clear tmp
end
load(fullfile(funPath,inDir,'featSel.mat'));



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
        error(['Please put this toolbox in Matlab path:' newline 'https://github.com/DrosteEffect/BrewerMap'])
    end
    save(filename,'cMap_F','cMap_vein');
end

f = {};
%% Plot Brain
f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
axBak = plotIm(axes,brain(:,:,slice));
title('BOLD image (1 TR)')


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
    end
end

if 1
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
            try
                featVals = featVal(:,featInd);
            catch
                keyboard
            end
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

if 0
    figure('WindowStyle','docked');
    corrplot(featVal(indIn,:),'varnames',featLabel)
end

%% Save figures
if 0
    filename = fullfile(funPath,outDir,mfilename);
    
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,[subjList{subjInd}]);
    for i = 1:length(f)
        curF = f{i};
        curF.Color = 'none';
        set(findobj(curF.Children,'type','Axes'),'color','none')
        curFile = [filename '_' num2str(i)];
        curExt = 'svg';
        saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
        curF.Color = 'w';
        curExt = 'fig';
        saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
        curExt = 'jpg';
        saveas(curF,[curFile '.' curExt]); if p.figOption.verbose; disp([curFile '.' curExt]); end
    end
end
