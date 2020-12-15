function defineAndShowMasks(fitType,threshType,veinPerc)
noMovement = 1;
actuallyRun = 1;
saveFig = 1;
plotAllSubj = 1;
if ~exist('veinPerc','var') || isempty(veinPerc)
    doVein = 1;
    veinSource = 'reducedModelResid'; % 'reducedModelResid' (stimulus-driven signal included in std) or 'fullModelResid (stimulus-driven signal excluded in std)'
    veinPerc = 20;
elseif veinPerc==0
    doVein = 0;
else
    doVein = 1;
    veinSource = 'reducedModelResid'; % 'reducedModelResid' (stimulus-driven signal included in std) or 'fullModelResid (stimulus-driven signal excluded in std)'
end
featSelContrast1 = 'anyCondActivation';
if ~exist('fitType','var') || isempty(fitType)
    fitType = 'fixed'; % 'mixed' (different regressors for each run) or 'fixed' (different regressors for each session)
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
funLevel3 = 'z';
if noMovement
    sinFitFile = 'v1SinCos_1perRun.mat';
    hrFitFile = 'v1resp_1perRun_resp.mat';
else
    sinFitFile = 'v1SinCos_1perRun_move12.mat';
    hrFitFile = 'v1resp_1perRun_move12_resp.mat';
end
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';

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
disp(['OUT: sinusoidal fit results (C-derived\DecodingHR\fun\' funLevel3 ')'])


for subjInd = 1:length(subjList)
    %% Load fun data
    load(fullfile(funPath,funLevel2,subjList{subjInd},sinFitFile),'results')
    tmp = load(fullfile(funPath,funLevel2,subjList{subjInd},hrFitFile),'results');
    resultsResp = tmp.results; clear tmp
    
    %% Define all labels
    runLabel = [results.OLS.mixed.inputs.opt.runLabel{:}];
    sessLabel = [results.OLS.mixed.inputs.opt.sessionLabel{:}];
    condLabel = repmat(1:3,[length(results.OLS.mixed.inputs.opt.runLabel)/3 1]); condLabel = condLabel(:)';
    
    % Repeat labels (tricky)
    %sort runLabels
    runLabelOrig = runLabel;
    [~,b] = sort(runLabel);
    runLabel = runLabel(b);
    sessLabel = sessLabel(b);
    condLabel = condLabel(b);
%     [runLabel' sessLabel' condLabel']
    %define repeat labels
    repeatLabel = nan(size(sessLabel));
    for sessInd = 1:2
        tmp = repmat(1:sum(sessLabel==sessInd)/3,3,1); tmp = tmp(:)';
        repeatLabel(sessLabel==sessInd) = tmp;
    end
%     [runLabel' sessLabel' condLabel' repeatLabel']
    %sort back to original order (important because later resorting may depend on that order)
    [~,b] = ismember(runLabelOrig,runLabel);
    runLabel = runLabel(b);
    sessLabel = sessLabel(b);
    condLabel = condLabel(b);
    repeatLabel = repeatLabel(b);
%     [runLabel' sessLabel' condLabel' repeatLabel']
        
    %% Set masks
	% Brain
    tmpFilename = dir(fullfile(funPath,funLevel1,subjList{subjInd},'trun101_preprocessed.nii*'));
    if all(size(tmpFilename)==[1 1]); tmpFilename = tmpFilename.name; else; error('X'); end
    a = load_nii(fullfile(funPath,funLevel1,subjList{subjInd},tmpFilename));
    brain = flip(permute(a.img,[3 1 2 4]),1); clear a
    brain = brain(:,:,:,1);

    % Mask of V1
    tmpFilename = dir(fullfile(anatPath,anatLevel,subjList{subjInd},'v1.nii*'));
    if all(size(tmpFilename)==[1 1]); tmpFilename = tmpFilename.name; else; error('X'); end
%     tmpFilename = ls(fullfile(anatPath,anatLevel,subjList{subjInd},'v1.nii*'));
    a = load_nii(fullfile(anatPath,anatLevel,subjList{subjInd},tmpFilename));
    maskV1 = flipdim(permute(a.img,[3 1 2 4]),1); clear a
    maskV1(:,:,1) = zeros(size(maskV1,1),size(maskV1,2));% Remove corrupted slices
    maskV1(:,:,end) = zeros(size(maskV1,1),size(maskV1,2));% Remove corrupted slices
    
    % Mask of eccentricity (stimulus retinotopic representation)
    tmpFilename = dir(fullfile(anatPath,anatLevel,subjList{subjInd},'lh.ecc.nii*'));
    if all(size(tmpFilename)==[1 1]); tmpFilename = tmpFilename.name; else; error('X'); end
    a = load_nii(fullfile(anatPath,anatLevel,subjList{subjInd},tmpFilename));
    eccL = flipdim(permute(a.img,[3 1 2 4]),1); clear a
    tmpFilename = dir(fullfile(anatPath,anatLevel,subjList{subjInd},'rh.ecc.nii*'));
    if all(size(tmpFilename)==[1 1]); tmpFilename = tmpFilename.name; else; error('X'); end
    a = load_nii(fullfile(anatPath,anatLevel,subjList{subjInd},tmpFilename));
    eccR = flipdim(permute(a.img,[3 1 2 4]),1); clear a
    find((eccL~=0 & eccR~=0));
    ecc = eccL; ecc(eccR~=0) = eccR(eccR~=0); clear eccL eccR
    maskECC = ecc>1 & ecc<6;
    
    % Mask of ROI
    maskROI = maskV1 & maskECC;
    
    % Mask of fit area
    maskFitArea = results.mask;
    
    % Mask of veins (session-specific)
    exclusion.subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
    exclusion.subj = 2;
    exclusion.sess = {1};
    exclusion.run = {5};
    exclusion.cond = {2};
    if doVein
        for sessInd = 1:2
            if length(exclusion.subj)>1; error('Need to code for multiple exclusions'); end
            % deal with exclusion here
            if strcmp(exclusion.subjList{exclusion.subj},subjList{subjInd}) && exclusion.sess{1}==sessInd
                ind = repeatLabel==exclusion.run{1} & sessLabel==exclusion.sess{1};
                ind = ~ind;
            else
                ind = true(size(runLabel));
            end
            vein.(['sess' num2str(sessInd)]).noiseOverMean = nan(size(maskROI));
            switch veinSource
                case 'reducedModelResid'
                    vein.(['sess' num2str(sessInd)]).noiseOverMean(maskFitArea) = mean(results.OLS.mixed.veinFull(:,:,:,ind),4);
                case 'fullModelResid'
                    vein.(['sess' num2str(sessInd)]).noiseOverMean(maskFitArea) = mean(results.OLS.mixed.veinReduced(:,:,:,ind),4);
            end
            vein.(['sess' num2str(sessInd)]).thresh = prctile(vein.(['sess' num2str(sessInd)]).noiseOverMean(maskROI),100-veinPerc);
            vein.(['sess' num2str(sessInd)]).mask = vein.(['sess' num2str(sessInd)]).noiseOverMean>vein.(['sess' num2str(sessInd)]).thresh;
        end
    end
    
    %% Set feature selection stats (session-specific)
    for sessInd = 1:2
        switch fitType
            case 'mixed'
                warning(['***Exclusion not performed on the random-model (mixed-model) stats' newline 'but should have only a minor impact on voxel selection for a signle session of a single subject***'])
                F = results.OLS.(fitType).(['Fsess' num2str(sessInd)]).val.F;
                df = results.OLS.(fitType).(['Fsess' num2str(sessInd)]).df;
            case 'fixed'
                F = results.OLS.(fitType).(['sess' num2str(sessInd)]).F.val.F;
                df = results.OLS.(fitType).(['sess' num2str(sessInd)]).F.df;
            otherwise
                error('X')
        end
        
        tmp = nan(size(brain)); tmp(maskFitArea) = F;
        F = tmp;
        [~,P] = getPfromF(F,df);
        [tmp,~] = getPfromF(F(maskROI),df);
        FDR = nan(size(F)); FDR(maskROI) = tmp;
        
        featSel.(['sess' num2str(sessInd)]).(featSelContrast1).F = F;
        featSel.(['sess' num2str(sessInd)]).(featSelContrast1).P = P;
        featSel.(['sess' num2str(sessInd)]).(featSelContrast1).FDR = FDR;
    end
    
    %% Plot masking
    if plotAllSubj || subjInd==1
        sess = 'sess1';
        slice = 11;
        
        
        %Get some nice colormaps
        filename = fullfile(pwd,mfilename);
        if ~exist(filename,'dir'); mkdir(filename); end
        filename = fullfile(filename,'cmap');
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
        
        f = [];
        %Brain
        f(end+1) = figure('WindowStyle','docked','color','w');
        axBak = plotIm(axes,brain(:,:,slice));
        title('BOLD image (1 TR)')
        
        %Activation overlay
        f(end+1) = figure('WindowStyle','docked','color','w');
        axBak = plotIm(axes,brain(:,:,slice));
        im = featSel.(sess).(featSelContrast1).F(:,:,slice);
        im = log(im);
        cLim = [min(im(:)) max(im(:))];
        cLim(1) = cLim(1) + diff(cLim)*0.2;
        axOver = plotIm(axes,im,cLim);
        alphaData = ~isnan(im);
        makeOverlay(axBak,axOver,alphaData,cMap_F,'log',cLim)
        ylabel('Activation Level (F value)');
        %thresholded
        f(end+1) = figure('WindowStyle','docked','color','w');
        axBak = plotIm(axes,brain(:,:,slice));
        maskFull = maskROI;
        if ~strcmp(threshType,'none')
            maskThresh = featSel.(sess).(featSelContrast1).(upper(threshType))<threshVal;
            maskFull = maskFull & maskThresh;
        end
        if doVein
            maskFull = maskFull & ~vein.(sess).mask;
        end
        im(~maskFull(:,:,slice)) = nan;
        axOver = plotIm(axes,im,cLim);
        alphaData = ~isnan(im);
        makeOverlay(axBak,axOver,alphaData,cMap_F,'log',cLim)
        ylabel({'Activation Level (F value)' ['\fontsize{8}Thresholded (' threshType '<' num2str(threshVal,'%0.2f') ') within ROI']});
        if ~strcmp(threshType,'none')
            %hist
            f(end+1) = figure('WindowStyle','docked','color','w');
            vol = featSel.(sess).(featSelContrast1).F(maskROI);
            vol = log(vol);
            hHist = histogram(vol,1000); hold on
            tip = 0.001;
            [~,bHigh] = min(abs(cumsum(hHist.Values)./sum(hHist.Values)-(1-tip)));
            [~,bLow] = min(abs(cumsum(hHist.Values)./sum(hHist.Values)-tip));
            xLim = [hHist.BinEdges(bLow) hHist.BinEdges(bHigh+1)];
            xlim(xLim);
            i = 0;
            done = 0;
            XTicks = [];
            XTickLabel = [];
            while ~done
                XTicks = [XTicks (0:0.1:1)*(10^i)];
                XTickLabel = [XTickLabel (10^i)];
                i = i+1;
                done = XTicks(end)>exp(xLim(2));
            end
            imMin = xLim(1);
            imMax = xLim(2);
            XTicks = sort(unique(XTicks));
            XTicks = XTicks(XTicks>exp(imMin) & XTicks<exp(imMax));
            ax = gca;
            ax.XTick = log(XTicks);
            ax.XTickLabel = cellstr(num2str(XTicks'));
            tmp = ax.XTickLabel;
            tmp(~ismember(XTicks,XTickLabel)) = {''};
            ax.XTickLabel = tmp;
            ax.XTickLabel(1) = {num2str(XTicks(1))};
            ax.XTickLabel(end) = {num2str(XTicks(end))};
            ax.TickDir = 'out';
            ax.Box = 'off';
            hHist.FaceColor = 'k'; hHist.EdgeColor = 'none';
            ylabel('voxel count')
            xlabel('F-value')
            tmp = featSel.(sess).(featSelContrast1).(upper(threshType))(maskROI);
            [~,b] = min(abs(tmp-threshVal));
            yLim = ylim;
            plot([1 1].*vol(b),yLim,'r')
            hTex = text(vol(b),yLim(2),['F=' num2str(exp(vol(b)),'%0.1f') ', ' threshType '=' num2str(threshVal,'%0.2f')]);
            hTex.VerticalAlignment = 'top';
            hTex.Position(1) = hTex.Position(1) + hTex.Extent(3).*0.02;
        end
        
        if doVein
            %Vein overlay
            f(end+1) = figure('WindowStyle','docked','color','w');
            axBak = plotIm(axes,brain(:,:,slice));
            axOver = plotIm(axes,vein.(sess).noiseOverMean(:,:,slice));
            alphaData = ~isnan(featSel.(sess).(featSelContrast1).F(:,:,slice));
            cLim = vein.(sess).noiseOverMean(:,:,slice);
            cLim = cLim(maskROI(:,:,slice));
            cLim = [min(cLim) max(cLim)];
            makeOverlay(axBak,axOver,alphaData,cMap_vein,'lin',cLim)
            ylabel('Veinness (signal std / signal mean)')
            %hist
            f(end+1) = figure('WindowStyle','docked','color','w');
            vol = vein.(sess).noiseOverMean(maskROI);
            hHist = histogram(vol,100); hold on
            tip = 0.001;
            [~,bHigh] = min(abs(cumsum(hHist.Values)./sum(hHist.Values)-(1-tip)));
            [~,bLow] = min(abs(cumsum(hHist.Values)./sum(hHist.Values)-tip));
            xLim = [hHist.BinEdges(bLow) hHist.BinEdges(bHigh+1)];
            xlim(xLim);
            ax = gca;
            ax.TickDir = 'out';
            ax.Box = 'off';
            hHist.FaceColor = 'k'; hHist.EdgeColor = 'none';
            ylabel('voxel count')
            xlabel('Veinness (BOLD std/mean)')
            tmp = vein.(sess).thresh;
            yLim = ylim;
            plot([1 1].*tmp,yLim,'b')
            hTex = text(tmp,yLim(2),num2str(100-veinPerc,'%0.0f%%ile'));
            hTex.VerticalAlignment = 'top';
            hTex.Position(1) = hTex.Position(1) + hTex.Extent(3).*0.02;
        end
        
        %V1 retinotopic representation + Vein
        f(end+1) = figure('WindowStyle','docked','color','w');
        axBak = plotIm(axes,brain(:,:,slice));
        axV1 = plotIm(axes,double(maskROI(:,:,slice)));
        makeOverlay(axBak,axV1,maskROI(:,:,slice),[0 0 0; 255 255 0]./255)
        if doVein
            tmp = double(vein.(sess).mask(:,:,slice));
            tmp(~maskROI(:,:,slice)) = 0;
            axVein = plotIm(axes,tmp);
            makeOverlay(axBak,axVein,tmp,[0 0 0; 0 140 225]./255)
        end
        
        % Save
        if saveFig
            filename = fullfile(pwd,mfilename);
            if ~exist(filename,'dir'); mkdir(filename); end
            filename = fullfile(filename,[subjList{subjInd}]);
            for i = 1:length(f)
                curF = figure(f(i));
                curF.Color = 'none';
                set(findobj(curF.Children,'type','Axes'),'color','none')
                saveas(curF,[filename '_' num2str(i) '.svg']); disp([filename '_' num2str(i) '.svg'])
                curF.Color = 'w';
%                 set(findobj(curF.Children,'type','Axes'),'color','w')
                saveas(curF,[filename '_' num2str(i)]); disp([filename '_' num2str(i) '.fig'])
                saveas(curF,[filename '_' num2str(i) '.jpg']); disp([filename '_' num2str(i) '.jpg'])
            end
        end
    end
    
    %% Stop here if not actually runing
    if ~actuallyRun
        return
    end
    
    %% Apply masks to fun data, vectorize voxels and split sessions and conditions
    % Split sessions and conditions
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
        ind(:) = maskROI(maskFitArea);
        
        % sin responses
        d.(sess).data = permute(data.(sess).data,[4 5 1 2 3]);
        d.(sess).data = d.(sess).data(:,:,ind);
        d.(sess).data = permute(d.(sess).data,[1 3 2]);
        d.(sess).hr = permute(hr.(sess),[4 5 6 1 2 3]);
        d.(sess).hr = d.(sess).hr(:,:,:,ind);
        d.(sess).hr = permute(d.(sess).hr,[2 4 3 1]);
        d.(sess).runLabel = permute(data.(sess).runLabel,[4 5 1 2 3]);
        d.(sess).runLabel = permute(d.(sess).runLabel,[1 3 2]);
        d.(sess).info = '%BOLD: run x vox x cond[ori1, ori2, plaid] X TR';
        
        % feature slection stats
        list = {'F' 'FDR' 'P'};
        for i = 1:3
            tmp = nan(sz(1:3));
            tmp(:) = featSel.(sess).(featSelContrast1).(list{i})(maskFitArea);
            featSelStats.(featSelContrast1).(sess).(list{i}) = tmp(ind)';
        end
        featSelStats.(featSelContrast1).(sess).info = '1 x vox';
        
        % vein
        if doVein
            list = {'noiseOverMean' 'mask'};
            for i = 1:length(list)
                tmp = nan(sz(1:3));
                tmp(:) = vein.(sess).(list{i})(maskFitArea);
                vein.(sess).(list{i}) = tmp(ind)';
            end
        else
            vein.(sess).mask = false(size(featSelStats.(featSelContrast1).(sess).F));
        end
        vein.(sess).info = '1 x vox';
        
    end
    clear hr
    
    %% Export some parameters
    param.subjList = subjList;
    param.featSelContrast1.name = featSelContrast1;
    param.featSelContrast1.threshType = threshType;
    param.featSelContrast1.threshVal = threshVal;
    if doVein
        param.vein.veinSource = veinSource;
        param.vein.veinPerc = veinPerc;
    end
    
    %% Save
    if ~exist(fullfile(funPath,funLevel3),'dir')
        mkdir(fullfile(funPath,funLevel3));
    end

    tmp = fullfile(funPath,funLevel3,[subjList{subjInd} '_' mfilename]);
    save(tmp,'d','featSelStats','vein','param')
    disp(['Saving to: ' tmp '.mat'])
end


