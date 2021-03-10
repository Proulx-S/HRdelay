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
featSelVal = nan([size(featSel2{subjInd,sessInd}.featVal,2) size(brain)]);
featSelValThresholded = nan([size(featSel2{subjInd,sessInd}.featVal,2) size(brain)]);
for i = 1:size(featSel2{subjInd,sessInd}.featVal,2)
    featSelVal(i,pAll{subjInd}.masks.roiMasks.v1) = featSel2{subjInd,sessInd}.featVal(:,i);
    featSelValThresholded(i,pAll{subjInd}.masks.roiMasks.v1) = featSel2{subjInd,sessInd}.featVal(:,i);
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

% Unthresholded maps
doIt = true(size(featSel2{subjInd,sessInd}.featLabel));
doIt(3) = false;
for featInd = 1:size(featSel2{subjInd,sessInd}.featLabel,2)
    if doIt(featInd)
        f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
        axBak = plotIm(axes,brain(:,:,slice));
        im = featSelVal(:,:,slice,featInd);
        switch featSel2{subjInd,sessInd}.featLabel{featInd}
            case 'veinScore'
                ylabel({'Vein score' 'BOLD var / BOLD mean'})
                cMap = cMap_vein;
                scale = 'log';
                im = log(im);
            case 'act'
                ylabel('Activation Level (F value)');
                cMap = cMap_F;
                scale = 'log';
                im = log(im);
            case 'discrim'
                ylabel('Discriminativeness (Hotelling''s T^2)');
                cMap = cMap_F;
                scale = 'lin';
            otherwise
                error('X')
        end
        cLim = [min(im(:)) max(im(:))];
%         cLim = prctile(im(:),[1 99]);
        cLim(1) = cLim(1) + diff(cLim)*0.2;
        axOver = plotIm(axes,im,cLim);
        alphaData = ~isnan(im);
        makeOverlay(axBak,axOver,alphaData,cMap,scale,cLim)
    end
end

% Thresholded maps
doIt = false(size(featSel2{subjInd,sessInd}.featLabel));
doIt(2) = true;
for featInd = 1:size(featSel2{subjInd,sessInd}.featLabel,2)
    if doIt(featInd)
        f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
        axBak = plotIm(axes,brain(:,:,slice));
        im = featSelValThresholded(:,:,slice,featInd);
        switch featSel2{subjInd,sessInd}.featLabel{featInd}
            case 'veinScore'
                ylabel({'Thresholded Vein score' 'BOLD var / BOLD mean'})
                cMap = cMap_vein;
                scale = 'log';
                im = log(im);
            case 'act'
                ylabel({'Thresholded Activation Level' '(F value)'});
                cMap = cMap_F;
                scale = 'log';
                im = log(im);
            case 'discrim'
                ylabel({'Thresholded Discriminativeness' '(Hotelling''s T^2)'});
                cMap = cMap_F;
                scale = 'lin';
            otherwise
                error('X')
        end
        cLim = [min(im(:)) max(im(:))];
%         cLim = prctile(im(:),[1 99]);
        cLim(1) = cLim(1) + diff(cLim)*0.2;
        axOver = plotIm(axes,im,cLim);
        alphaData = ~isnan(im);
        makeOverlay(axBak,axOver,alphaData,cMap,scale,cLim)
    end
end


% Histograms
doIt = true(size(featSel2{subjInd,sessInd}.featLabel));
for featInd = 1:size(featSel2{subjInd,sessInd}.featLabel,2)
    if doIt(featInd)
        f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
        tmp = log(featSelVal(:,:,:,featInd));
        [N,edges] = histcounts(tmp(v1));
        Nout = histcounts(tmp(v1 & ~featSelIn),edges);
        Nin = histcounts(tmp(v1 & featSelIn),edges);
        ctrs = (edges(1:end-1)+edges(2:end))/2; % Calculate the bin centers
        bar(ctrs, [Nout' Nin'],'stacked','FaceAlpha',0.6); hold on
        legend({'exlcuded voxels' 'included voxels'},'box','off');
        ax = gca;
        ax.Box = 'off';
%         delta = diff(exp(edges([1 end])));
%         delta = 1/10 .^ (-floor(log10(abs(delta))));
%         ax.XTick = log(floor(exp(edges(1))./delta):delta/10:ceil(exp(edges(1))./delta));
%         ax.XTickLabel = exp(ax.XTick);
%         XTickLabel = cellstr(ax.XTickLabel);
%         XTickLabel(~ismember(ax.XTick,log(floor(exp(edges(1))./delta):delta:ceil(exp(edges(1))./delta)))) = {''};
%         ax.XTickLabel = XTickLabel;
        ax.TickDir = 'out';
        switch featSel2{subjInd,sessInd}.featLabel{featInd}
            case 'veinScore'
                xlabel({'Vein score' 'BOLD var / BOLD mean'})
            case 'act'
                xlabel({'Activation Level' 'log(F value)'});
%                 max(featSel2{subjInd,sessInd}.featP{featInd}(featSel2{subjInd,sessInd}.indIn))
                [~,b] = min(abs(featSel2{subjInd,sessInd}.featP{featInd}-0.05));
                tmpF = log(featSel2{subjInd,sessInd}.featVal(b,featInd));
                plot([1 1].*tmpF,ylim,'k')
            case 'discrim'
                xlabel({'Discriminativeness' '(Hotelling''s T^2)'});
            otherwise
                error('X')
        end
    end
end

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
% featInd = ismember(featSel2{subjInd,sessInd}.featLabel,'act');
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
% featInd = ismember(featSel2{subjInd,sessInd}.featLabel,'veinScore');
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
% x = featSel2{subjInd,sessInd}.featVal(:,1);
% y = featSel2{subjInd,sessInd}.featVal(:,2);
% z = featSel2{subjInd,sessInd}.featVal(:,3);
% indIn = featSel2{subjInd,sessInd}.indIn;
% featLabel = featSel2{subjInd,sessInd}.featLabel;
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






