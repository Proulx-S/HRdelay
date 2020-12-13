function maskSinAndHrFit(fitType,threshType)
noMovement = 1;
actuallyRun = 1;
saveFig = 0;
plotAllSubj = 0;
doVein = 1;
veinPerc = 10;
if ~exist('fitType','var') || isempty(fitType)
    fitType = 'mixed'; % 'mixed' (different regressors for each run) or 'fixed' (different regressors for each session)
end
if ~exist('threshType','var') || isempty(threshType) % for plotting purpose only (does not affect data that is saved)
    threshType = 'p'; % 'none', 'p' or 'fdr'
end
threshVal = 0.05;
warning('off','images:imshow:magnificationMustBeFitForDockedFigure');

%% Define paths
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
dataDir = 'C-derived\DecodingHR';
anatPath = fullfile(repoPath,dataDir,'anat');
anatLevel = 'z';
funPath = fullfile(repoPath,dataDir,'fun');
funLevel1 = 'x';
funLevel2 = 'y';
funLevel3Sin = 'zSin';
funLevel3Hr = 'zHr';
if noMovement
    sinFitFile = 'v1SinCos_1perRun.mat';
    hrFitFile = 'v1resp_1perRun_resp.mat';
else
    sinFitFile = 'v1SinCos_1perRun_move12.mat';
    hrFitFile = 'v1resp_1perRun_move12_resp.mat';
end
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';
% subjList = {'02jp' '03sk' '04sp'}';

%make sure everything is forward slash for mac, linux pc compatibility
for tmpPath = {'repoPath' 'dataDir' 'anatPath' 'funPath'}
    eval([char(tmpPath) '(strfind(' char(tmpPath) ',''\''))=''/'';']);
end


if ~actuallyRun
    disp('Not actually running to save some time')
end
disp(['IN: anatomical V1 roi (C-derived\DecodingHR\anat\' anatLevel ')'])
disp(['IN: voxel visual field eccentricity (C-derived\DecodingHR\anat\' anatLevel ')'])
disp(['IN: sinusoidal fit results (C-derived\DecodingHR\fun\' funLevel2 ')'])
disp('F(IN)=OUT: masks the fit according to voxel eccentricity and activation level')
disp('Figures are additionally thresholded for activation level, but not the data that is saved!')
disp(['OUT: sinusoidal fit results (C-derived\DecodingHR\fun\' funLevel3Sin ')'])


for subjInd = 1:length(subjList)
    %% Load fun data
    load(fullfile(funPath,funLevel2,subjList{subjInd},sinFitFile),'results')
    tmp = load(fullfile(funPath,funLevel2,subjList{subjInd},hrFitFile),'results');
    resultsResp = tmp.results; clear tmp
    
    % brain
    tmpFilename = dir(fullfile(funPath,funLevel1,subjList{subjInd},'trun101_preprocessed.nii*'));
    if all(size(tmpFilename)==[1 1]); tmpFilename = tmpFilename.name; else; error('X'); end
    a = load_nii(fullfile(funPath,funLevel1,subjList{subjInd},tmpFilename));
    brain = flip(permute(a.img,[3 1 2 4]),1); clear a
    brain = brain(:,:,:,1);
    
    %% Set masks
    % maskV1
    tmpFilename = dir(fullfile(anatPath,anatLevel,subjList{subjInd},'v1.nii*'));
    if all(size(tmpFilename)==[1 1]); tmpFilename = tmpFilename.name; else; error('X'); end
%     tmpFilename = ls(fullfile(anatPath,anatLevel,subjList{subjInd},'v1.nii*'));
    a = load_nii(fullfile(anatPath,anatLevel,subjList{subjInd},tmpFilename));
    maskV1 = flipdim(permute(a.img,[3 1 2 4]),1); clear a
    maskV1(:,:,1) = zeros(size(maskV1,1),size(maskV1,2));% Remove corrupted slices
    maskV1(:,:,end) = zeros(size(maskV1,1),size(maskV1,2));% Remove corrupted slices
    
    % maskECC
    tmpFilename = dir(fullfile(anatPath,anatLevel,subjList{subjInd},'lh.ecc.nii*'));
    if all(size(tmpFilename)==[1 1]); tmpFilename = tmpFilename.name; else; error('X'); end
%     tmpFilename = ls(fullfile(anatPath,anatLevel,subjList{subjInd},'lh.ecc.nii*'));
    a = load_nii(fullfile(anatPath,anatLevel,subjList{subjInd},tmpFilename));
    eccL = flipdim(permute(a.img,[3 1 2 4]),1); clear a
    tmpFilename = dir(fullfile(anatPath,anatLevel,subjList{subjInd},'rh.ecc.nii*'));
    if all(size(tmpFilename)==[1 1]); tmpFilename = tmpFilename.name; else; error('X'); end
%     tmpFilename = ls(fullfile(anatPath,anatLevel,subjList{subjInd},'rh.ecc.nii*'));
    a = load_nii(fullfile(anatPath,anatLevel,subjList{subjInd},tmpFilename));
    eccR = flipdim(permute(a.img,[3 1 2 4]),1); clear a
    find((eccL~=0 & eccR~=0));
    ecc = eccL; ecc(eccR~=0) = eccR(eccR~=0); clear eccL eccR
    maskECC = ecc>1 & ecc<6;
    
    % maskAnat
    maskAnat = maskV1 & maskECC;
    
    % maskFit area
    maskFit = results.mask;
    
    % maskVein
    sessLabel = [results.inputs.opt.sessionLabel{:}]';
    for sessInd = 1:2
        vein.(['sess' num2str(sessInd)]).noiseOverMean = nan(size(maskAnat));
        vein.(['sess' num2str(sessInd)]).noiseOverMean(maskFit) = mean(results.OLS.mixed.vein(:,:,:,sessLabel==sessInd),4);
        vein.(['sess' num2str(sessInd)]).thresh = prctile(vein.(['sess' num2str(sessInd)]).noiseOverMean(maskAnat),100-veinPerc);
        vein.(['sess' num2str(sessInd)]).mask = vein.(['sess' num2str(sessInd)]).noiseOverMean>vein.(['sess' num2str(sessInd)]).thresh;
    end
    
    % feature selection
    for sessInd = 1:2
        switch fitType
            case 'mixed'
                eval(['F = results.OLS.mixed.Fsess' num2str(sessInd) '.val.F;']);
%                 eval(['F = results.OLS.mixed.Fsess' num2str(sessInd) '.val.F(maskAnat(maskFit));']);
                eval(['df = results.OLS.mixed.Fsess' num2str(sessInd) '.df;']);
            case 'fixed'
                eval(['F = results.OLS.fixed.sess' num2str(sessInd) '.F.val.F;']);
%                 eval(['F = results.OLS.fixed.sess' num2str(sessInd) '.F.val.F(maskAnat(maskFit));']);
                eval(['df = results.OLS.fixed.sess' num2str(sessInd) '.F.df;']);
            otherwise
                error('X')
        end
        
        tmp = nan(size(brain)); tmp(maskFit) = F;
        F = tmp;
        [~,P] = getPfromF(F,df);
        [tmp,~] = getPfromF(F(maskAnat),df);
        FDR = nan(size(F)); FDR(maskAnat) = tmp;
        
        featSel.(['sess' num2str(sessInd)]).anyCondActivation.F = F;
        featSel.(['sess' num2str(sessInd)]).anyCondActivation.P = P;
        featSel.(['sess' num2str(sessInd)]).anyCondActivation.FDR = FDR;
    end
    
    %% Plot masking
    if plotAllSubj || subjInd==1
        sess = 'sess1';
        slice = 11;
        
        
        %Get some noce colormaps
        filename = fullfile(pwd,mfilename);
        if ~exist(filename,'dir'); mkdir(filename); end
        filename = fullfile(filename,'cmap');
%         FMap = brewermap(256,'reds');
%         VeinMap = brewermap(256,'blues');
%         save(filename,'FMap','VeinMap');
        load(filename,'FMap','VeinMap');
        
        %Brain
        f = figure('WindowStyle','docked','color','w');
        axBak = plotIm(axes,brain(:,:,slice));
        title('BOLD image (1 TR)')
        
        %Activation overlay
        f = figure('WindowStyle','docked','color','w');
        axBak = plotIm(axes,brain(:,:,slice));
        axOver = plotIm(axes,log(featSel.(sess).anyCondActivation.F(:,:,slice)));
        alphaData = ~isnan(featSel.(sess).anyCondActivation.F(:,:,slice));
        makeOverlay(axOver,alphaData,FMap,'log')

        %Vein overlay
        f = figure('WindowStyle','docked','color','w');
        axBak = plotIm(axes,brain(:,:,slice));
        axOver = plotIm(axes,vein.(sess).noiseOverMean(:,:,slice));
        alphaData = ~isnan(featSel.(sess).anyCondActivation.F(:,:,slice));
        cLim = vein.(sess).noiseOverMean(:,:,slice);
        cLim = cLim(maskAnat(:,:,slice));
        cLim = [min(cLim) max(cLim)];
        makeOverlay(axOver,alphaData,VeinMap,'lin',cLim)
        
        
        %V1 retinotopic representation + Vein
        f = figure('WindowStyle','docked','color','w');
        axBak = plotIm(axes,brain(:,:,slice));
        axV1 = plotIm(axes,double(maskAnat(:,:,slice)));
        makeOverlay(axV1,maskAnat(:,:,slice),[0 0 0; 255 255 0]./255)
        tmp = double(vein.(sess).mask(:,:,slice));
        tmp(~maskAnat(:,:,slice)) = 0;
        axVein = plotIm(axes,tmp);
        makeOverlay(axVein,tmp,[0 0 0; 0 140 225]./255)
        
        
        
        
%         figure('WindowStyle','docked');
%         im = brain(:,:,slice);
%         im = uint8(im./max(im(:)).*255);
%         im = uint8(ind2rgb(im,gray(256)).*255);
%         im = permute(im,[3 1 2]);
%         
%         im(1,maskAnat(:,:,slice)) = 255;
%         im(2,maskAnat(:,:,slice)) = 255;
%         im(3,maskAnat(:,:,slice)) = 40;
%         
%         im(1,vein.(sess).mask(:,:,slice)) = 28;
%         im(2,vein.(sess).mask(:,:,slice)) = 150;
%         im(3,vein.(sess).mask(:,:,slice)) = 255;
%         
%         im = permute(im,[2 3 1]);
%         imshow(im)
%         title({'ROI + Veins' '\fontsize{8} (V1, 1 to 6 dva eccentricity)'})
%         
% %         ax3 = subplot(2,2,3); axis off
%         figure('WindowStyle','docked');
%         F = featSel.(sess).anyCondActivation.F(:,:,slice);
%         F(isnan(F)) = 0;
%         switch threshType
%             case 'p'
%                 thresh = featSel.(sess).anyCondActivation.P(:,:,slice);
%             case 'fdr'
%                 thresh = featSel.(sess).anyCondActivation.FDR(:,:,slice);
%             case 'none'
%             otherwise
%                 error('x')
%         end
%         if ~strcmp(threshType,'none')
%             F(thresh>threshVal) = 0; % functional threshold here
%         end
%         tmp = unique(F(:));
%         [~,b] = sort(tmp);
%         Fmin =  tmp(b(2)); clear b
%         Fmax = max(F(:));
%         imF = uint8((F-Fmin)./(Fmax-Fmin)*255);
%         imF = uint8(ind2rgb(imF,autumn(256))*255);
%         
%         im = brain(:,:,slice);
%         im = uint8(im./max(im(:)).*255);
%         im = uint8(ind2rgb(im,gray(256)).*255);
%         
%         im = permute(im,[3 1 2]);
%         imF = permute(imF,[3 1 2]);
%         for rgb = 1:3
%             im(rgb,F~=0) = imF(rgb,F~=0);
%         end
%         im = permute(im,[2 3 1]); clear imF
%         imshow(im)
%         
%         
%         
%         
%         
%         
%         
%         %vein overlay
%         f = figure('WindowStyle','docked','color','w');
%         axBak = plotIm(axes,brain(:,:,slice));
%         axOver = plotIm(axes,vein.(sess).noiseOverMean(:,:,slice));
%         alphaData = ~isnan(featSel.(sess).anyCondActivation.F(:,:,slice));
%         makeOverlay(axOver,alphaData,FMap,'log')
% 
%         ax = axF;
%         hIm = findobj(ax.Children,'Type','Image');
%         i = 0;
%         done = 0;
%         YTicks = [];
%         YTickLabel = [];
%         while ~done
%             YTicks = [YTicks (0:0.1:1)*(10^i)];
%             YTickLabel = [YTickLabel (10^i)];
%             i = i+1;
%             done = YTicks(end)>exp(axF.CLim(2));
%         end
%         imMin = min(hIm.CData(:));
%         imMax = max(hIm.CData(:));
%         YTicks = sort(unique(YTicks));
%         YTicks = YTicks(YTicks>exp(imMin) & YTicks<exp(imMax));
%         ax.YTick = interp1([imMin imMax],[1 size(im,1)],log(YTicks));
%         ax.YTickLabel = cellstr(num2str(YTicks'));
%         tmp = ax.YTickLabel;
%         tmp(~ismember(YTicks,YTickLabel)) = {''};
%         ax.YTickLabel = tmp;
%         ax.YTickLabel(1) = {num2str(YTicks(1))};
%         ax.YTickLabel(end) = {num2str(YTicks(end))};
%         ax.TickDir = 'out';
%         
%         
%         ylabel('F value')
%         title({'Visual Activation' '\fontsize{8} (within ROI)'})
%         
%         
%         
%         
%         
%         
%         
%         axF.Children.AlphaData(logical(axF.Children.AlphaData)) = ~isnan(featSel.(sess).anyCondActivation.F(:,:,slice));
%         axF.Children.AlphaData = axF.Children.AlphaData.*1
%         
%         
%         
%         im = featSel.(sess).anyCondActivation.F(:,:,slice);
%         cLim = [min(im(:)) max(im(:))];
%         cb = linspace(cLim(1),cLim(2),size(im,1))';
%         cb = repmat(cb,1,round(size(im,1).*cbWidth));
%         im = cat(2,cb,im);
%         axF.Children.CData = im;
%         axF.CLim = cLim;
%         axF.Colormap = FMap;
%         
%         hIm.YData = fliplr(hIm.YData);
%         axBrain = gca;
%         axBrain.YDir = 'normal';
%         axBrain.PlotBoxAspectRatio = [1+cbWidth 1 1];
%         xticks([]); yticks([]);
%         axBrain.Colormap = colormap('gray');
%         
%         
%         
%         
%         
%         
%         
%         F = featSel.(sess).anyCondActivation.F(:,:,slice);
%         alphaData = ~isnan(F);
%         F = log(F);
%         Fmin = min(F(:));
%         Fmax = max(F(:));
%         imF = uint8((F-Fmin)./(Fmax-Fmin)*255);
%         imF = uint8(ind2rgb(imF,FMap)*255);
%         axF.Children.CData = imF;
%         axF.Children.AlphaData = alphaData.*1;
%         axF.PositionConstraint = 'innerposition';
%         cb=colorbar(axF)
%         
%         
%         
%         tmp = colormap('gray')
%         axBrain = gca;
%         axBrain.PlotBoxAspectRatio = [1 1 1];
%         
%         im = uint8(im./max(im(:)).*255);
%         im = uint8(ind2rgb(im,gray(256)).*255);
%         
%         colormap(size(im,1))
%         %add colorbar
%         imCB = uint8(autumn(size(im,1))*255);
%         imCB = permute(flipud(imCB),[1 3 2]);
%         imCB = repmat(imCB,[1 round(size(im,1).*0.03) 1]);
%         im = cat(2,imCB,im);
%         hIm = imshow(im);
%         
%         
%         imshow(im);
%         axBrain = gca;
%         
%         axF = copyobj(axBrain,f);
%         F = featSel.(sess).anyCondActivation.F(:,:,slice);
%         alphaData = ~isnan(F);
%         F = log(F);
%         Fmin = min(F(:));
%         Fmax = max(F(:));
%         imF = uint8((F-Fmin)./(Fmax-Fmin)*255);
%         imF = uint8(ind2rgb(imF,FMap)*255);
%         axF.Children.CData = imF;
%         axF.Children.AlphaData = alphaData.*1;
%         axF.PositionConstraint = 'innerposition';
%         cb=colorbar(axF)
%         cb
%         
%         
%         
%         
%         
%         
%         
%             fSubj(subjInd).Color = 'none';
%             set(findobj(fSubj(subjInd).Children,'type','Axes'),'color','none')
%             saveas(fSubj(subjInd),[filename '.svg']); disp([filename '.svg'])
%             fSubj(subjInd).Color = 'w';
%             set(findobj(fSubj(subjInd).Children,'type','Axes'),'color','w')
%             saveas(fSubj(subjInd),filename); disp([filename '.fig'])
%         save()
%         
%         
%         %Activation
%         figure('WindowStyle','docked');
%         F = log(featSel.(sess).anyCondActivation.F(:,:,slice));
%         Fmin = min(F(:));
%         Fmax = max(F(:));
%         imF = uint8((F-Fmin)./(Fmax-Fmin)*255);
%         imF = uint8(ind2rgb(imF,autumn(256))*255);
%         
%         im = brain(:,:,slice);
%         im = uint8(im./max(im(:)).*255);
%         im = uint8(ind2rgb(im,gray(256)).*255);
%         
%         im = permute(im,[3 1 2]);
%         imF = permute(imF,[3 1 2]);
%         for rgb = 1:3
%             im(rgb,~isnan(F)) = imF(rgb,~isnan(F));
%         end
%         im = permute(im,[2 3 1]); clear imF
%         %add colorbar
%         imCB = uint8(autumn(size(im,1))*255);
%         imCB = permute(flipud(imCB),[1 3 2]);
%         imCB = repmat(imCB,[1 round(size(im,1).*0.03) 1]);
%         im = cat(2,imCB,im);
%         hIm = imshow(im);
%         
%         ax = gca;
%         hIm.YData = [128 1];
%         ax.YDir = 'normal';
%         ax.YAxis.Visible = 'on';
%         YTicks = sort(unique([0:0.1:1 1:1:10 10:10:100 100:100:1000]));
%         YTicks = YTicks(YTicks>exp(Fmin) & YTicks<exp(Fmax));
%         ax.YTick = interp1([Fmin Fmax],[1 size(im,1)],log(YTicks));
%         ax.YTickLabel = cellstr(num2str(YTicks'));
%         YTickLabel = ax.YTickLabel;
%         YTickLabel(~ismember(YTicks,[1 10 100 1000])) = {''};
%         ax.YTickLabel = YTickLabel;
%         ax.YTickLabel(1) = {num2str(YTicks(1))};
%         ax.YTickLabel(end) = {num2str(YTicks(end))};
%         ylabel('F value')
%         title({'Visual Activation' '\fontsize{8} (within ROI)'})
%         
%         
%         %V1 + Vein
% %         ax2 = subplot(2,2,2);
%         figure('WindowStyle','docked');
%         im = brain(:,:,slice);
%         im = uint8(im./max(im(:)).*255);
%         im = uint8(ind2rgb(im,gray(256)).*255);
%         im = permute(im,[3 1 2]);
%         
%         im(1,maskAnat(:,:,slice)) = 255;
%         im(2,maskAnat(:,:,slice)) = 255;
%         im(3,maskAnat(:,:,slice)) = 40;
%         
%         im(1,vein.(sess).mask(:,:,slice)) = 28;
%         im(2,vein.(sess).mask(:,:,slice)) = 150;
%         im(3,vein.(sess).mask(:,:,slice)) = 255;
%         
%         im = permute(im,[2 3 1]);
%         imshow(im)
%         title({'ROI + Veins' '\fontsize{8} (V1, 1 to 6 dva eccentricity)'})
%         
% %         ax3 = subplot(2,2,3); axis off
%         figure('WindowStyle','docked');
%         F = featSel.(sess).anyCondActivation.F(:,:,slice);
%         F(isnan(F)) = 0;
%         switch threshType
%             case 'p'
%                 thresh = featSel.(sess).anyCondActivation.P(:,:,slice);
%             case 'fdr'
%                 thresh = featSel.(sess).anyCondActivation.FDR(:,:,slice);
%             case 'none'
%             otherwise
%                 error('x')
%         end
%         if ~strcmp(threshType,'none')
%             F(thresh>threshVal) = 0; % functional threshold here
%         end
%         tmp = unique(F(:));
%         [~,b] = sort(tmp);
%         Fmin =  tmp(b(2)); clear b
%         Fmax = max(F(:));
%         imF = uint8((F-Fmin)./(Fmax-Fmin)*255);
%         imF = uint8(ind2rgb(imF,autumn(256))*255);
%         
%         im = brain(:,:,slice);
%         im = uint8(im./max(im(:)).*255);
%         im = uint8(ind2rgb(im,gray(256)).*255);
%         
%         im = permute(im,[3 1 2]);
%         imF = permute(imF,[3 1 2]);
%         for rgb = 1:3
%             im(rgb,F~=0) = imF(rgb,F~=0);
%         end
%         im = permute(im,[2 3 1]); clear imF
%         imshow(im)
% 
%         axCB = copyobj(ax3,gcf);
%         axes(axCB);
%         CB = permute(flipud(autumn(size(ax1.Children.CData,1))),[1 3 2]);
%         CB = repmat(CB,[1 round(0.1*size(ax1.Children.CData,2)) 1]);
%         imshow(CB);
%         axCB.Position(1) = ax3.Position(1) - ax3.Position(3)*0.5;
%         axCB.YAxis.Visible = 'on';
%         axCB.YAxis.TickValues = [1 size(ax1.Children.CData,1)];
%         if strcmp(threshType,'none')
%             axCB.YAxis.Label.String = 'F';
%             axCB.YAxis.TickLabels = {num2str(Fmax,'%0.1f') num2str(Fmin,'%0.1f')};
%         else
%             tmp = threshVal-thresh(:); tmp(tmp<0) = nan; [~,b] = min(tmp);
%             axCB.YAxis.TickLabels = {num2str(max(F(:)),'%0.1f') num2str(F(b),'%0.1f')};
%             if strcmp(threshType,'p')
%                 axCB.YAxis.Label.String = {'F' ['p<' num2str(threshVal,'%0.2f')]};
%             elseif strcmp(threshType,'fdr')
%                 axCB.YAxis.Label.String = {'F' ['FDR<' num2str(threshVal,'%0.2f')]};
%             end
%         end
%         axCB.YAxis.Label.Rotation = 90;
%         axCB.YAxis.Label.VerticalAlignment = 'top';
%         axCB.YAxis.Label.HorizontalAlignment = 'center';
%         title('Activation Map')
%         
%         
%         
%         figure('WindowStyle','docked');
%         ax4 = subplot(2,2,4);
%         hm = histogram(featSel.(sess).anyCondActivation.F(maskAnat),1000); hold on
%         xlim([0 20])
%         ylim(ylim*1.05)
%         hm.FaceColor = 'k';
%         ax4.PlotBoxAspectRatio = [1 1 1];
%         ylabel('vox count');
%         xlabel('F')
%         switch threshType
%             case 'none'
%             case 'p'
%                 line([1 1].*F(b),ylim,'Color','r');
%                 text(F(b)+diff(xlim)*0.01,ax4.YLim(2)*0.95,['p=' num2str(threshVal,'%0.2f')])
%             case 'fdr'
%                 line([1 1].*F(b),ylim,'Color','r');
%                 text(F(b)+diff(xlim)*0.01,ax4.YLim(2)*0.95,['FDR=' num2str(threshVal,'%0.2f')])
%             otherwise
%                 error('X')
%         end
%         title('All ROI voxels')
%         if ~strcmp(threshType,'none')
%             ax4.XTick = sort([ax4.XTick round(F(b),2)]);
%         end
%         ax4.TickDir = 'out';
%         ax4.XTick(ax4.XTick==0) = [];
        
        
        
        
        if saveFig
            filename = fullfile(pwd,mfilename);
            if ~exist(filename,'dir'); mkdir(filename); end
            filename = fullfile(filename,[subjList{subjInd}]);
            fSubj(subjInd).Color = 'none';
            set(findobj(fSubj(subjInd).Children,'type','Axes'),'color','none')
            saveas(fSubj(subjInd),[filename '.svg']); disp([filename '.svg'])
            fSubj(subjInd).Color = 'w';
            set(findobj(fSubj(subjInd).Children,'type','Axes'),'color','w')
            saveas(fSubj(subjInd),filename); disp([filename '.fig'])
        end
    end
    
    %% Stop here if not actually runing
    if ~actuallyRun
        return
    end
    
    %% Apply masks to fun data, vectorize voxels and split sessions and conditions
    % Split sessions and conditions
    runLabel = [results.inputs.opt.runLabel{:}];
    sessLabel = [results.inputs.opt.sessionLabel{:}];
    condLabel = repmat(1:3,[length(results.inputs.opt.runLabel)/3 1]); condLabel = condLabel(:)';
    [X,Y] = pol2cart(results.OLS.mixed.delay,results.OLS.mixed.amp);
    for sessInd = 1:2
        sess = ['sess' num2str(sessInd)];
        ind = sessLabel==sessInd;
        data.(sess).data = complex(X(:,:,:,ind),Y(:,:,:,ind));
        data.(sess).condLabel = permute(condLabel(ind),[1 3 4 2]);
        data.(sess).runLabel = permute(runLabel(ind),[1 3 4 2]);
        
        cond1 = data.(sess).condLabel==1;
        cond2 = data.(sess).condLabel==2;
        cond3 = data.(sess).condLabel==3;
        data.(sess).data = cat(5,data.(sess).data(:,:,:,cond1),data.(sess).data(:,:,:,cond2),data.(sess).data(:,:,:,cond3));
        data.(sess).runLabel = cat(5,data.(sess).runLabel(:,:,:,cond1),data.(sess).runLabel(:,:,:,cond2),data.(sess).runLabel(:,:,:,cond3));
        data.(sess).info = 'x X y X z X run X cond[ori1,ori2,plaid]';
        
        data.(sess) = rmfield(data.(sess),'condLabel');
    end
    clear X Y
    
    hr.sess1 = resultsResp.OLS.mixed.sess1.resp;
    hr.sess2 = resultsResp.OLS.mixed.sess2.resp;
    hr.info = 'x X y X z X TR X run X cond[ori1,ori2,plaid]';
    
    % Mask anat and vectorize
    for sessInd = 1:2
        sess = ['sess' num2str(sessInd)];
        sz = size(data.(sess).data);
        ind = false(sz(1:3));
        ind(:) = maskAnat(maskFit);
        
        % sin responses
        d.(sess).data = permute(data.(sess).data,[4 5 1 2 3]);
        d.(sess).data = d.(sess).data(:,:,ind);
        d.(sess).data = permute(d.(sess).data,[1 3 2]);
        d.(sess).runLabel = permute(data.(sess).runLabel,[4 5 1 2 3]);
        d.(sess).runLabel = permute(d.(sess).runLabel,[1 3 2]);
        d.(sess).info = '%BOLD: run x vox x cond[ori1, ori2, plaid]';
        
        % stats
        list = {'F' 'FDR' 'P'};
        for i = 1:3
            tmp = nan(sz(1:3));
            tmp(:) = featSel.(sess).anyCondActivation.(list{i})(maskFit);
            stats.(sess).(list{i}) = tmp(ind)';
        end
        stats.(sess).info = '1 x vox';
        
        % response shape
        hrTmp.(sess).data = permute(hr.(sess),[4 5 6 1 2 3]);
        hrTmp.(sess).data = hrTmp.(sess).data(:,:,:,ind);
        hrTmp.(sess).data = permute(hrTmp.(sess).data,[2 4 3 1]);
        hrTmp.(sess).info = '%BOLD: run X vox X cond[ori1,ori2,plaid] X TR';
        
        % vein
        list = {'noiseOverMean' 'mask'};
        for i = 1:length(list)
            tmp = nan(sz(1:3));
            tmp(:) = vein.(sess).(list{i})(maskFit);
            vein.(sess).(list{i}) = tmp(ind)';
        end
        vein.(sess).info = '1 x vox';
    end
    hr = hrTmp; clear hrTmp
    
    
    %% Save
    if ~exist(fullfile(funPath,funLevel3Sin),'dir')
        mkdir(fullfile(funPath,funLevel3Sin));
    end
    if ~exist(fullfile(funPath,funLevel3Hr),'dir')
        mkdir(fullfile(funPath,funLevel3Hr));
    end
    
    if noMovement
        save(fullfile(funPath,funLevel3Sin,[subjList{subjInd} '_' mfilename '_noMovement']),'d','stats','vein')
        save(fullfile(funPath,funLevel3Hr,[subjList{subjInd} '_' mfilename '_noMovement']),'hr','stats','vein')
    else
        save(fullfile(funPath,funLevel3Sin,[subjList{subjInd} '_' mfilename]),'d','stats','vein')
        save(fullfile(funPath,funLevel3Hr,[subjList{subjInd} '_' mfilename]),'hr','stats','vein')
    end
end


