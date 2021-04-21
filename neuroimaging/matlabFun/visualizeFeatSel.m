function visualizeFeatSel(p)
if ~isfield(p,'figOption') || isempty(p.figOption)
    p.figOption.verbose = 1;
    p.figOption.subjInd = 1;
    p.figOption.sessInd = 1;
end


%% Define paths
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'd';
            inDir2  = 'b';
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
%     if subjInd==1
%         tmp = load(fullfile(funPath,inDir2,[subjList{subjInd} '.mat']),'p','d');
%     else
        tmp = load(fullfile(funPath,inDirX,[subjList{subjInd} '.mat']),'p');
%     end
    pAll{subjInd} = tmp.p; clear tmp
end
load(fullfile(funPath,inDir,'featSel.mat'));



% pAll{1}
% featSel{1}
% featSel{1}

%% Initialize subj, sess and masks
subjInd = 1;
sessInd = 1;
for i = 1:length(featSel{subjInd,sessInd}.condPairList)
    if length(featSel{subjInd,sessInd}.condPairList{i})==3 ...
            && all(featSel{subjInd,sessInd}.condPairList{i} == [1 2 3])
        condPairInd = i;
    end
end
brain = pAll{subjInd}.brain.mean(:,:,:,sessInd);
roi = pAll{subjInd}.voxProp.area.ind==pAll{subjInd}.voxProp.area.indList(ismember(pAll{subjInd}.voxProp.area.labelList,p.featSel.fov.areaLabel))...
    & ~pAll{subjInd}.voxProp.censorMask;


featSelIn = false(size(brain));
featSelIn(roi) = featSel{subjInd,sessInd}.indIn(:,:,condPairInd);
featSelSeqIn = false([length(featSel{subjInd,sessInd}.featSeq.featSelList) size(brain)]);
featSelSeqIn(:,roi) = featSel{subjInd,sessInd}.featSeq.featIndIn(:,:,condPairInd)';
featSelSeqIn = permute(featSelSeqIn,[2 3 4 1]);
featSelVal = nan([size(featSel{subjInd,sessInd}.featSeq.featVal,2) size(brain)]);
featSelValThresholded = nan([size(featSel{subjInd,sessInd}.featSeq.featVal,2) size(brain)]);
featInfo = featSel{subjInd,sessInd}.featSeq.featSelList;
condInfo = featSel{subjInd,sessInd}.featSeq.condPairList(:,:,condPairInd);
featLabel = cell(size(featInfo));
for featInd = 1:size(featInfo,2)
    tmp = strsplit(featInfo{featInd},': ');
    featLabel(featInd) = tmp(1);
    featSelVal(featInd,roi) = featSel{subjInd,sessInd}.featSeq.featVal(:,featInd,condPairInd);
    featSelValThresholded(featInd,roi) = featSel{subjInd,sessInd}.featSeq.featVal(:,featInd,condPairInd);
    featSelValThresholded(featInd,~featSelIn) = nan;
end
featSelVal = permute(featSelVal,[2 3 4 1]);
featSelValThresholded = permute(featSelValThresholded,[2 3 4 1]);





sess = ['sess' num2str(sessInd)];
slice = 11;
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
% doIt(ismember(featLabel,'retinoFov')) = false;
% doIt(3) = false;

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

%% Plot response vector distribution
tmp = load(fullfile(funPath,'d',[subjList{subjInd} '.mat']));
d = tmp.res.(sess);
f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
x = mean(d.sin(:,:,:,:),4);
tmp = [];
featInd = ismember(featLabel,'respVecSig');
% indPrev = featSel{subjInd,sessInd}.featSeq.featPrevInd(:,featInd,condIndPairInd);
ind = featSel{subjInd,sessInd}.featSeq.featIndIn(:,featInd,condPairInd);
polarscatter(angle(x(ind)),abs(x(ind)),eps,'.k'); hold on
polarscatter(angle(x(~ind)),abs(x(~ind)),eps,'.r');

% f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
% x = wrapToPi(angle(x)-angle(mean(x)))+angle(mean(x));
% x = -x./pi.*6;
% figure('WindowStyle','docked');
% histogram(x,100)
% ind
% xlabel('delay (sec)')
% xlim([0 12])



%% Plot response vector maps
% tmp = load(fullfile(funPath,'c',[subjList{subjInd} '.mat']));
% masks = tmp.p.masks; clear tmp
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
%     cLim = [-1 1].*max(abs(im(:)));
    cLim = [cLim(1) + diff(cLim)*0.2 cLim(2) - diff(cLim)*0.2];
    axOver = plotIm(axes,im,cLim);
    alphaData = ~isnan(im);
    makeOverlay(axBak,axOver,alphaData,cMap,scale,cLim)
    ylabel([yLabel {[num2str(cLim(1)) ' to ' num2str(cLim(2))]}])
% end


%% Plot histograms
doIt = true(size(featInfo));
% doIt(ismember(featLabel,'retinoFov')) = false;
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
                %                 xlabel({'Vein score' 'BOLD var / BOLD mean'})    
            case 'vein'
                xLabel = {'Vein score' 'log( BOLD var / BOLD mean )'};
                qtile = 1-thresh.percentile/100;
                scale = 'log';
                xLim = 'auto';
%                 xlabel({'Vein score' 'BOLD var / BOLD mean'})
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
        
%         [nnz(tmpInd(roi) & ~featSelIn(roi))
%         nnz(~tmpInd(roi) & featSelIn(roi))
%         nnz(tmpInd(roi) & featSelIn(roi))
%         nnz(~tmpInd(roi) & ~featSelIn(roi))]
%         nnz(roi)
%         
%         indOut_thisFeat = ~featSelSeqIn(:,:,:,featInd);
%         indOut_otherFeat = ~featSelSeqIn; indOut_otherFeat(:,:,:,featInd) = []; indOut_otherFeat = any(indOut_otherFeat,4);
%         indIn = featSelIn;
%         
%         [nnz(indOut_thisFeat(roi) & ~indOut_otherFeat(roi))
%         nnz(~indOut_thisFeat(roi) & indOut_otherFeat(roi))
%         nnz(indOut_thisFeat(roi) & indOut_otherFeat(roi))
%         nnz(~indOut_thisFeat(roi) & ~indOut_otherFeat(roi))]
%         nnz(roi(roi))
%         nnz(indIn(roi))
%         
%         [~,edges] = histcounts(tmp);
%         N1 = histcounts(tmp(roi & indOut_thisFeat & ~indOut_otherFeat),edges)';
%         N2 = histcounts(tmp(roi & ~indOut_thisFeat & indOut_otherFeat),edges)';
%         N3 = histcounts(tmp(roi & indOut_thisFeat & indOut_otherFeat),edges)';
%         N4 = histcounts(tmp(roi & ~indOut_thisFeat & ~indOut_otherFeat),edges)';
%         ctrs = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
%         hBar = bar(ctrs, [N4 N1 N2 N3],1,'stacked','FaceAlpha',0.6); hold on
%         legend({'4' '1' '2' '3'})
        
        
        
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
            featVals = featSel{subjInd,sessInd}.featSeq.featVal(:,featInd,condPairInd);
            qtiles = featSel{subjInd,sessInd}.featSeq.featQtile(:,featInd,condPairInd);
            ps = featSel{subjInd,sessInd}.featSeq.featP(:,featInd,condPairInd);
            fdrs = nan(size(ps));
            startInds = featSel{subjInd,sessInd}.featSeq.featIndStart(:,featInd,condPairInd);
            switch thresh.threshMethod
                case '%ile'
                    featVal = interp1(qtiles(~isnan(qtiles)&~isnan(featVals)),featVals(~isnan(qtiles)&~isnan(featVals)),qtile);
                case 'p'
                    featVal = interp1(ps,featVals,thresh.threshVal);
                case 'fdr'
                    fdrs(startInds) = mafdr(ps(startInds),'BHFDR',true);
                    [~,b] = sort(abs(fdrs-thresh.threshVal));
                    featVal = mean(featVals(b(1:2)));
                otherwise
                    error('X')
            end
            switch scale
                case 'log'
                    featVal = log(featVal);
                case 'lin'
            end
            hThresh = plot([1 1].*featVal,ylim,':k');
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
        
        
        if 0
            switch featLabel{featInd}
                case 'vein'
                    switch thresh.threshMethod
                        case '%ile'
                            featVals = featSel{subjInd,sessInd}.featSeq.featVal(:,featInd);
                            qtiles = featSel{subjInd,sessInd}.featSeq.featQtile(:,featInd);
                            featVal = interp1(qtiles,featVals,qtile);
                            switch scale
                                case 'log'
                                    featVal = log(featVal);
                                case 'lin'
                            end
                            plot([1 1].*featVal,ylim,'r');
                        case 'p'
                            error('code that')
                        case 'fdr'
                            error('code that')
                        otherwise
                            error('X')
                    end
                case 'act'
                    xlabel({'Activation Level' 'log( F value )'});
                    %                 max(featSel{subjInd,sessInd}.featP{featInd}(featSel{subjInd,sessInd}.indIn(:,:,condPairInd)))
                    [~,b] = min(abs(featSel{subjInd,sessInd}.featP{featInd}-0.05));
                    tmpF = log(featSel{subjInd,sessInd}.featVal(b,featInd));
                    plot([1 1].*tmpF,ylim,'k')
                case 'respVecSig'
                    xlabel({'Response Vector Strenght' 'log( Hotelling''s U )'});
                    %                 max(featSel{subjInd,sessInd}.featP{featInd}(featSel{subjInd,sessInd}.indIn(:,:,condPairInd)))
                    [~,b] = min(abs(featSel{subjInd,sessInd}.featP{featInd}-0.05));
                    tmpF = log(featSel{subjInd,sessInd}.featVal(b,featInd));
                    plot([1 1].*tmpF,ylim,'k')
                case 'respVecDiff'
                    xlabel({'Discriminativeness' 'log( Hotelling''s U )'});
                otherwise
                    error('X')
            end
        end
    end
end

figure('WindowStyle','docked');
featVal = log(featSel{subjInd,sessInd}.featSeq.featVal(featSel{subjInd,sessInd}.indIn(:,:,condPairInd),:,condPairInd));
corrplot(featVal,'varnames',featLabel)
% figure('WindowStyle','docked');
% featVal = log(featSel{subjInd,sessInd}.featSeq.featVal(:,:,condPairInd));
% corrplot(featVal,'varnames',featLabel)


if 0
    % Save
    if p.figOption.save
        filename = fullfile(pwd,mfilename);
        if ~exist(filename,'dir'); mkdir(filename); end
        filename = fullfile(filename,[subjList{subjInd}]);
        for i = 1:length(f)
            curF = f{i};
            curF.Color = 'none';
            set(findobj(curF.Children,'type','Axes'),'color','none')
            curFile = [filename '_' num2str(i)];
            curExt = 'svg';
            saveas(curF,[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
            curF.Color = 'w';
            %                 set(findobj(curF.Children,'type','Axes'),'color','w')
            curExt = 'fig';
            saveas(curF,[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
            curExt = 'jpg';
            saveas(curF,[curFile '.' curExt]); if verbose; disp([curFile '.' curExt]); end
        end
    end
end





% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % Activation overlay
% f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
% axBak = plotIm(axes,brain(:,:,slice));
% featInd = ismember(featInfo,'act');
% im = featSelVal(:,:,slice,featInd);
% % im = featSel.(sess).anyCondActivation_F(:,:,slice);
% im = log(im);
% cLim = [min(im(:)) max(im(:))];
% cLim(1) = cLim(1) + diff(cLim)*0.2;
% axOver = plotIm(axes,im,cLim);
% alphaData = ~isnan(im);
% makeOverlay(axBak,axOver,alphaData,cMap_F,'log',cLim)
% ylabel('Activation Level (F value)');
% % thresholded
% f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
% axBak = plotIm(axes,brain(:,:,slice));
% im = featSelValThresholded(:,:,slice,2);
% % im = featSel.(sess).anyCondActivation_F(:,:,slice);
% im = log(im);
% cLim = [min(im(:)) max(im(:))];
% cLim(1) = cLim(1) + diff(cLim)*0.2;
% axOver = plotIm(axes,im,cLim);
% alphaData = ~isnan(im);
% makeOverlay(axBak,axOver,alphaData,cMap_F,'log',cLim)
% ylabel('Activation Level (thresholded F value)');
% 
% % % hist
% % f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
% % vol = featSel.(sess).anyCondActivation_F(maskROI);
% % vol = log(vol);
% % hHist = histogram(vol,'BinWidth',0.1); hold on
% % tip = 0.001;
% % [~,bHigh] = min(abs(cumsum(hHist.Values)./sum(hHist.Values)-(1-tip)));
% % [~,bLow] = min(abs(cumsum(hHist.Values)./sum(hHist.Values)-tip));
% % xLim = [hHist.BinEdges(bLow) hHist.BinEdges(bHigh+1)];
% % xlim(xLim);
% % i = 0;
% % done = 0;
% % XTicks = [];
% % XTickLabel = [];
% % while ~done
% %     XTicks = [XTicks (0:0.1:1)*(10^i)];
% %     XTickLabel = [XTickLabel (10^i)];
% %     i = i+1;
% %     done = XTicks(end)>exp(xLim(2));
% % end
% % imMin = xLim(1);
% % imMax = xLim(2);
% % XTicks = sort(unique(XTicks));
% % XTicks = XTicks(XTicks>exp(imMin) & XTicks<exp(imMax));
% % ax = gca;
% % ax.XTick = log(XTicks);
% % ax.XTickLabel = cellstr(num2str(XTicks'));
% % tmp = ax.XTickLabel;
% % tmp(~ismember(XTicks,XTickLabel)) = {''};
% % ax.XTickLabel = tmp;
% % ax.XTickLabel(1) = {num2str(XTicks(1))};
% % ax.XTickLabel(end) = {num2str(XTicks(end))};
% % ax.TickDir = 'out';
% % ax.Box = 'off';
% % hHist.FaceColor = 'k'; hHist.EdgeColor = 'none';
% % ylabel('voxel count')
% % xlabel('F-value')
% % drawnow
% % if ~strcmp(featSel.(sess).anyCondActivation_threshType,'none')
% %     tmp = featSel.(sess).(['anyCondActivation_' (upper(featSel.(sess).anyCondActivation_threshType))])(maskROI);
% %     [~,b] = min(abs(tmp-featSel.(sess).anyCondActivation_threshVal));
% %     yLim = ylim;
% %     plot([1 1].*vol(b),yLim,'r')
% %     drawnow
% %     hTex = text(vol(b),yLim(2),['F=' num2str(exp(vol(b)),'%0.1f') ', ' featSel.(sess).anyCondActivation_threshType '=' num2str(featSel.(sess).anyCondActivation_threshVal,'%0.2f')]);
% %     hTex.VerticalAlignment = 'top';
% %     hTex.Position(1) = hTex.Position(1) + hTex.Extent(3).*0.02;
% %     drawnow
% % end
% 
% % Vein overlay
% f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
% axBak = plotIm(axes,brain(:,:,slice));
% featInd = ismember(featInfo,'veinScore');
% im = log(featSelVal(:,:,slice,featInd));
% axOver = plotIm(axes,im);
% alphaData = v1(:,:,slice);
% cLim = [min(im(:)) max(im(:))];
% makeOverlay(axBak,axOver,alphaData,cMap_vein,'log',cLim)
% % axOver = plotIm(axes,featSel.(sess).vein_score(:,:,slice));
% % alphaData = ~isnan(featSel.(sess).vein_score(:,:,slice));
% % cLim = featSel.(sess).vein_score(:,:,slice);
% % cLim = cLim(maskROI(:,:,slice));
% % cLim = [min(cLim) max(cLim)];
% % makeOverlay(axBak,axOver,alphaData,cMap_vein,'lin',cLim)
% ylabel({'Veinness' 'signal std / signal mean'})
% % %hist
% % f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
% % vol = featSel.(sess).vein_score(maskROI);
% % hHist = histogram(vol,'BinWidth',0.002); hold on
% % tip = 0.001;
% % [~,bHigh] = min(abs(cumsum(hHist.Values)./sum(hHist.Values)-(1-tip)));
% % [~,bLow] = min(abs(cumsum(hHist.Values)./sum(hHist.Values)-tip));
% % xLim = [hHist.BinEdges(bLow) hHist.BinEdges(bHigh+1)];
% % xlim(xLim);
% % ax = gca;
% % ax.TickDir = 'out';
% % ax.Box = 'off';
% % hHist.FaceColor = 'k'; hHist.EdgeColor = 'none';
% % ylabel('voxel count')
% % xlabel('Vein score (BOLD std/mean)')
% % drawnow
% % if featSel.(sess).vein_doIt
% %     tmp = featSel.(sess).vein_thresh;
% %     yLim = ylim;
% %     plot([1 1].*tmp,yLim,'b')
% %     hTex = text(tmp,yLim(2),num2str(100-veinPerc,'%0.0f%%ile'));
% %     hTex.VerticalAlignment = 'top';
% %     hTex.Position(1) = hTex.Position(1) + hTex.Extent(3).*0.02;
% % end
% 
% % %V1 retinotopic representation + Vein
% % f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
% % axBak = plotIm(axes,brain(:,:,slice));
% % axV1 = plotIm(axes,double(maskROI(:,:,slice)));
% % makeOverlay(axBak,axV1,maskROI(:,:,slice),[0 0 0; 255 255 0]./255)
% % if featSel.(sess).vein_doIt
% %     tmp = double(featSel.(sess).vein_mask(:,:,slice));
% %     tmp(~maskROI(:,:,slice)) = 0;
% %     axVein = plotIm(axes,tmp);
% %     makeOverlay(axBak,axVein,tmp,[0 0 0; 0 140 225]./255)
% % end
% 
% 
% x = featSel{subjInd,sessInd}.featVal(:,1);
% y = featSel{subjInd,sessInd}.featVal(:,2);
% z = featSel{subjInd,sessInd}.featVal(:,3);
% indIn = featSel{subjInd,sessInd}.indIn(:,:,condPairInd);
% featLabel = featInfo;
% 
% figure('WindowStyle','docked');
% scatter3(x(indIn),y(indIn),z(indIn),'k.'); hold on
% scatter3(x(~indIn),y(~indIn),z(~indIn),'r.');
% % ax = gca;
% % ax.CameraPosition = camPos;
% xlabel(['' featLabel{1} ''])
% ylabel(['' featLabel{2} ''])
% zlabel(['' featLabel{3} ''])
% ax = gca;
% ax.XAxis.Scale = 'log';
% ax.YAxis.Scale = 'log';
% ax.ZAxis.Scale = 'log';
% %     xlabel(['zscore(-log(' info1{1+1} '))'])
% %     ylabel(['zscore(log(' info1{1+2} '))'])
% %     zlabel(['zscore(log(' info1{1+3} '))'])
% legend({'include' 'exclude'})






