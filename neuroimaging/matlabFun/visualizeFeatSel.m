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
    if subjInd==1
        tmp = load(fullfile(funPath,inDir2,[subjList{subjInd} '.mat']),'p');
    else
        tmp = load(fullfile(funPath,inDirX,[subjList{subjInd} '.mat']),'p');
    end
    pAll{subjInd} = tmp.p;
end
load(fullfile(funPath,inDir,'featSel.mat'));



pAll
featSel{1}
featSel2{1}

subjInd = 1;
sessInd = 1;
brain = pAll{subjInd}.brain.mean(:,:,:,sessInd);
v1 = pAll{subjInd}.masks.roiMasks.v1;
featSelIn = false(size(brain));
featSelIn(pAll{subjInd}.masks.roiMasks.v1) = featSel2{subjInd,sessInd}.indIn;
featSelVal = nan([size(featSel2{subjInd,sessInd}.percAll,2) size(brain)]);
featSelValThresholded = nan([size(featSel2{subjInd,sessInd}.percAll,2) size(brain)]);
for i = 1:size(featSel2{subjInd,sessInd}.percAll,2)
    featSelVal(i,pAll{subjInd}.masks.roiMasks.v1) = featSel2{subjInd,sessInd}.percAll(:,i);
    featSelValThresholded(i,pAll{subjInd}.masks.roiMasks.v1) = featSel2{subjInd,sessInd}.percAll(:,i);
    featSelValThresholded(i,~featSelIn) = nan;
end
featSelVal = permute(featSelVal,[2 3 4 1]);
featSelValThresholded = permute(featSelValThresholded,[2 3 4 1]);





sess = ['sess' num2str(sessInd)];
slice = 11;
visibility = true;

% Get some nice colormaps
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
% Brain
f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
axBak = plotIm(axes,brain(:,:,slice));
title('BOLD image (1 TR)')

% Activation overlay
f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
axBak = plotIm(axes,brain(:,:,slice));
featInd = ismember(featSel2{subjInd,sessInd}.featLabel,'act');
im = featSelVal(:,:,slice,featInd);
% im = featSel.(sess).anyCondActivation_F(:,:,slice);
im = log(im);
cLim = [min(im(:)) max(im(:))];
cLim(1) = cLim(1) + diff(cLim)*0.2;
axOver = plotIm(axes,im,cLim);
alphaData = ~isnan(im);
makeOverlay(axBak,axOver,alphaData,cMap_F,'log',cLim)
ylabel('Activation Level (F value)');
% thresholded
f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
axBak = plotIm(axes,brain(:,:,slice));
im = featSelValThresholded(:,:,slice,2);
% im = featSel.(sess).anyCondActivation_F(:,:,slice);
im = log(im);
cLim = [min(im(:)) max(im(:))];
cLim(1) = cLim(1) + diff(cLim)*0.2;
axOver = plotIm(axes,im,cLim);
alphaData = ~isnan(im);
makeOverlay(axBak,axOver,alphaData,cMap_F,'log',cLim)
ylabel('Activation Level (thresholded F value)');

% % hist
% f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
% vol = featSel.(sess).anyCondActivation_F(maskROI);
% vol = log(vol);
% hHist = histogram(vol,'BinWidth',0.1); hold on
% tip = 0.001;
% [~,bHigh] = min(abs(cumsum(hHist.Values)./sum(hHist.Values)-(1-tip)));
% [~,bLow] = min(abs(cumsum(hHist.Values)./sum(hHist.Values)-tip));
% xLim = [hHist.BinEdges(bLow) hHist.BinEdges(bHigh+1)];
% xlim(xLim);
% i = 0;
% done = 0;
% XTicks = [];
% XTickLabel = [];
% while ~done
%     XTicks = [XTicks (0:0.1:1)*(10^i)];
%     XTickLabel = [XTickLabel (10^i)];
%     i = i+1;
%     done = XTicks(end)>exp(xLim(2));
% end
% imMin = xLim(1);
% imMax = xLim(2);
% XTicks = sort(unique(XTicks));
% XTicks = XTicks(XTicks>exp(imMin) & XTicks<exp(imMax));
% ax = gca;
% ax.XTick = log(XTicks);
% ax.XTickLabel = cellstr(num2str(XTicks'));
% tmp = ax.XTickLabel;
% tmp(~ismember(XTicks,XTickLabel)) = {''};
% ax.XTickLabel = tmp;
% ax.XTickLabel(1) = {num2str(XTicks(1))};
% ax.XTickLabel(end) = {num2str(XTicks(end))};
% ax.TickDir = 'out';
% ax.Box = 'off';
% hHist.FaceColor = 'k'; hHist.EdgeColor = 'none';
% ylabel('voxel count')
% xlabel('F-value')
% drawnow
% if ~strcmp(featSel.(sess).anyCondActivation_threshType,'none')
%     tmp = featSel.(sess).(['anyCondActivation_' (upper(featSel.(sess).anyCondActivation_threshType))])(maskROI);
%     [~,b] = min(abs(tmp-featSel.(sess).anyCondActivation_threshVal));
%     yLim = ylim;
%     plot([1 1].*vol(b),yLim,'r')
%     drawnow
%     hTex = text(vol(b),yLim(2),['F=' num2str(exp(vol(b)),'%0.1f') ', ' featSel.(sess).anyCondActivation_threshType '=' num2str(featSel.(sess).anyCondActivation_threshVal,'%0.2f')]);
%     hTex.VerticalAlignment = 'top';
%     hTex.Position(1) = hTex.Position(1) + hTex.Extent(3).*0.02;
%     drawnow
% end

% Vein overlay
f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
axBak = plotIm(axes,brain(:,:,slice));
featInd = ismember(featSel2{subjInd,sessInd}.featLabel,'veinScore');
im = log(featSelVal(:,:,slice,featInd));
axOver = plotIm(axes,im);
alphaData = v1(:,:,slice);
cLim = [min(im(:)) max(im(:))];
makeOverlay(axBak,axOver,alphaData,cMap_vein,'log',cLim)
% axOver = plotIm(axes,featSel.(sess).vein_score(:,:,slice));
% alphaData = ~isnan(featSel.(sess).vein_score(:,:,slice));
% cLim = featSel.(sess).vein_score(:,:,slice);
% cLim = cLim(maskROI(:,:,slice));
% cLim = [min(cLim) max(cLim)];
% makeOverlay(axBak,axOver,alphaData,cMap_vein,'lin',cLim)
ylabel({'Veinness' 'signal std / signal mean'})
% %hist
% f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
% vol = featSel.(sess).vein_score(maskROI);
% hHist = histogram(vol,'BinWidth',0.002); hold on
% tip = 0.001;
% [~,bHigh] = min(abs(cumsum(hHist.Values)./sum(hHist.Values)-(1-tip)));
% [~,bLow] = min(abs(cumsum(hHist.Values)./sum(hHist.Values)-tip));
% xLim = [hHist.BinEdges(bLow) hHist.BinEdges(bHigh+1)];
% xlim(xLim);
% ax = gca;
% ax.TickDir = 'out';
% ax.Box = 'off';
% hHist.FaceColor = 'k'; hHist.EdgeColor = 'none';
% ylabel('voxel count')
% xlabel('Vein score (BOLD std/mean)')
% drawnow
% if featSel.(sess).vein_doIt
%     tmp = featSel.(sess).vein_thresh;
%     yLim = ylim;
%     plot([1 1].*tmp,yLim,'b')
%     hTex = text(tmp,yLim(2),num2str(100-veinPerc,'%0.0f%%ile'));
%     hTex.VerticalAlignment = 'top';
%     hTex.Position(1) = hTex.Position(1) + hTex.Extent(3).*0.02;
% end

% %V1 retinotopic representation + Vein
% f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
% axBak = plotIm(axes,brain(:,:,slice));
% axV1 = plotIm(axes,double(maskROI(:,:,slice)));
% makeOverlay(axBak,axV1,maskROI(:,:,slice),[0 0 0; 255 255 0]./255)
% if featSel.(sess).vein_doIt
%     tmp = double(featSel.(sess).vein_mask(:,:,slice));
%     tmp(~maskROI(:,:,slice)) = 0;
%     axVein = plotIm(axes,tmp);
%     makeOverlay(axBak,axVein,tmp,[0 0 0; 0 140 225]./255)
% end

% Save
if 0
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







