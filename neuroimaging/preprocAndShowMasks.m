function preprocAndShowMasks(featSel_bSess,figOption,verbose)
close all
if ~exist('verbose','var')
    verbose = 1;
end
noMovement = 1;
if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end
% Between-session feature selection
if ~exist('featSel_bSess','var') || isempty(featSel_bSess)
    featSel_bSess = struct('activation',[],'vein',[],'discrim',[]);
end
% activation
if ~isfield(featSel_bSess,'activation') || isempty(featSel_bSess.activation)
    featSel_bSess.activation.doIt = 0;
end
if featSel_bSess.activation.doIt
    if ~isfield(featSel_bSess.activation,'fitType') || isempty(featSel_bSess.activation.fitType)
        featSel_bSess.activation.fitType = 'fixed';
    end
    if ~isfield(featSel_bSess.activation,'threshType') || isempty(featSel_bSess.activation.threshType)
        featSel_bSess.activation.threshType = 'p';
    end
end
% vein
if ~isfield(featSel_bSess,'vein') || isempty(featSel_bSess.vein)
    featSel_bSess.vein.doIt = 0;
end
if featSel_bSess.vein.doIt
    if ~isfield(featSel_bSess.vein,'percentile') || isempty(featSel_bSess.vein.percentile)
        featSel_bSess.vein.percentile = 20;
    end
    if ~isfield(featSel_bSess.vein,'fullModelResid') || isempty(featSel_bSess.vein.fullModelResid)
        featSel_bSess.vein.fullModelResid = 20;
    end
end
% discrim
if ~isfield(featSel_bSess,'discrim') || isempty(featSel_bSess.discrim)
    featSel_bSess.discrim.doIt = 0;
end
if featSel_bSess.discrim.doIt
    if ~isfield(featSel_bSess.discrim,'percentile') || isempty(featSel_bSess.discrim.percentile)
        featSel_bSess.vein.percentile = 20;
    end
end


doVein = featSel_bSess.vein.doIt;
veinSource = featSel_bSess.vein.source; % 'reducedModelResid' (stimulus-driven signal included in std) or 'fullModelResid (stimulus-driven signal excluded in std)'
veinPerc = featSel_bSess.vein.percentile;

fitType = featSel_bSess.activation.fitType; % 'mixed' (different regressors for each run) or 'fixed' (different regressors for each session)
threshType = featSel_bSess.activation.threshType; % 'none', 'p' or 'fdr'
threshVal = featSel_bSess.activation.threshVal;
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

    %% Brain image
    tmpFilename = dir(fullfile(funPath,funLevel1,subjList{subjInd},'trun101_preprocessed.nii*'));
    if all(size(tmpFilename)==[1 1]); tmpFilename = tmpFilename.name; else; error('X'); end
    a = load_nii(fullfile(funPath,funLevel1,subjList{subjInd},tmpFilename));
    brain = flip(permute(a.img,[3 1 2 4]),1); clear a
    brain = brain(:,:,:,1);
%     figure('WindowStyle','docked')
%     imagesc(brain(:,:,10)); colormap gray
    %% ROI mask
    % Get mask of V1
    tmpFilename = dir(fullfile(anatPath,anatLevel,subjList{subjInd},'v1.nii*'));
    if all(size(tmpFilename)==[1 1]); tmpFilename = tmpFilename.name; else; error('X'); end
    %     tmpFilename = ls(fullfile(anatPath,anatLevel,subjList{subjInd},'v1.nii*'));
    a = load_nii(fullfile(anatPath,anatLevel,subjList{subjInd},tmpFilename));
    maskV1 = flipdim(permute(a.img,[3 1 2 4]),1); clear a
    maskV1(:,:,1) = zeros(size(maskV1,1),size(maskV1,2));% Remove corrupted slices
    maskV1(:,:,end) = zeros(size(maskV1,1),size(maskV1,2));% Remove corrupted slices
%     figure('WindowStyle','docked')
%     imagesc(maskV1(:,:,10)); colormap gray

    % Get mask of eccentricity (stimulus retinotopic representation)
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
%     figure('WindowStyle','docked')
%     imagesc(maskECC(:,:,10)); colormap gray

    % Combine the two
    maskROI = maskV1 & maskECC;
%     figure('WindowStyle','docked')
%     imagesc(maskROI(:,:,10)); colormap gray

    % Also get th mask of fit area (because not all voxels were fitted)
    maskFitArea = results.mask;
%     figure('WindowStyle','docked')
%     imagesc(maskFitArea(:,:,10)); colormap gray
    

    %% Split sessions and conditions
    [X,Y] = pol2cart(results.OLS.mixed.delay,results.OLS.mixed.amp);
    for sessInd = 1:2
        sess = ['sess' num2str(sessInd)];
        runInd = sessLabel==sessInd;
        data.(sess).data = complex(X(:,:,:,runInd),Y(:,:,:,runInd));
        data.(sess).condLabel = permute(condLabel(runInd),[1 3 4 2]);
        data.(sess).runLabel = permute(runLabel(runInd),[1 3 4 2]);
        %         [squeeze(data.(sess).condLabel) squeeze(data.(sess).runLabel)]
        
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
    hr.info = 'x X y X z X TR X repeat X cond[ori1,ori2,plaid]';
    
    %     figure('WindowStyle','docked')
    %     imagesc(abs(data.(sess).data(:,:,10,1,1)))
    
    
    %% Stats for later between-session feature selection
    exclusion.subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
    exclusion.subj = 2;
    exclusion.sess = {1};
    exclusion.run = {4};
    exclusion.cond = {2};
    for sessInd = 1:2
        sess = ['sess' num2str(sessInd)];
        % First deal with exclusion
        if length(exclusion.subj)>1; error('Need to code for multiple exclusions'); end
        if strcmp(exclusion.subjList{exclusion.subj},subjList{subjInd}) && exclusion.sess{1}==sessInd

            runInd = repeatLabel==exclusion.run{1} & sessLabel==exclusion.sess{1};
            runInd = ~runInd;
        else
            runInd = true(size(runLabel));
        end
        
        
        % Activated voxels
        if ~featSel_bSess.activation.doIt
            error('code that')
        else
            % get activation map
            featSel.(sess).anyCondActivation_doIt = ~strcmp(threshType,'none');
            featSel.(sess).anyCondActivation_fitType = fitType;
            switch fitType
                case 'mixed'
                    warning(['***Exclusion not performed on the random-model (mixed-model) stats' newline 'but should have only a minor impact on voxel selection for a signle session of a single subject***'])
                    F = results.OLS.(fitType).(['Fsess' num2str(sessInd)]).val.F;
                    df = results.OLS.(fitType).(['Fsess' num2str(sessInd)]).df;
                case 'fixed'
                    F = results.OLS.(fitType).(sess).F.val.F;
                    df = results.OLS.(fitType).(sess).F.df;
                otherwise
                    error('X')
            end
            tmp = nan(size(brain)); tmp(maskFitArea) = F;
            F = tmp;
            [~,P] = getPfromF(F,df);
            [tmp,~] = getPfromF(F(maskROI),df);
            FDR = nan(size(F)); FDR(maskROI) = tmp;
            featSel.(sess).anyCondActivation_F = F;
            featSel.(sess).anyCondActivation_P = P;
            featSel.(sess).anyCondActivation_FDR = FDR;
            featSel.(sess).anyCondActivation_threshType = threshType;
            featSel.(sess).anyCondActivation_threshVal = threshVal;
            % get activation mask
            if featSel.(sess).anyCondActivation_doIt
                featSel.(sess).anyCondActivation_mask = featSel.(sess).(['anyCondActivation_' upper(featSel.(sess).anyCondActivation_threshType)])<featSel.(sess).anyCondActivation_threshVal;
            else
                featSel.(sess).anyCondActivation_mask = true(size(featSel.(sess).anyCondActivation_F));
            end
%             figure('WindowStyle','docked')
%             imagesc(featSel.(sess).anyCondActivation_F(:,:,10))
%             imagesc(featSel.(sess).anyCondActivation_mask(:,:,10))
%             imagesc(maskROI(:,:,10))
%             imagesc(featSel.(sess).anyCondActivation_mask(:,:,10)&maskROI(:,:,10))
        end
        
        % Vein voxels
        featSel.(sess).vein_doIt = doVein;
        featSel.(sess).vein_source = veinSource;
        featSel.(sess).vein_score = nan(size(brain));
        featSel.(sess).vein_perc = veinPerc;
        featSel.(sess).vein_thresh = nan;
        featSel.(sess).vein_mask = false(size(brain));
        if ~featSel_bSess.vein.doIt
            error('code that')
        else
            % get vein score
            switch veinSource
                case 'fullModelResid'
                    featSel.(sess).vein_score(maskFitArea) = mean(results.OLS.mixed.veinFull(:,:,:,runInd & sessLabel==sessInd),4);
                case 'reducedModelResid'
                    featSel.(sess).vein_score(maskFitArea) = mean(results.OLS.mixed.veinReduced(:,:,:,runInd & sessLabel==sessInd),4);
            end
            % get vein threshold as the percentile of ACTIVE (see Activation map above) voxels in ROI
            maskTmp = featSel.(sess).anyCondActivation_mask & maskROI;
            %             % get vein threshold as percentile of voxels in ROI
            %             maskTmp = maskROI;
            featSel.(sess).vein_thresh = prctile(featSel.(sess).vein_score(maskTmp),100-featSel.(sess).vein_perc);
            % get vein mask
            featSel.(sess).vein_mask = featSel.(sess).vein_score>featSel.(sess).vein_thresh;
%             figure('WindowStyle','docked')
%             im = featSel.(sess).vein_score;
%             imagesc(im(:,:,10))
%             maskTmp = ~featSel.(sess).vein_mask;
%             imagesc(maskTmp(:,:,10));
%             maskTmp = ~featSel.(sess).vein_mask & maskROI;
%             imagesc(maskTmp(:,:,10));
%             maskTmp = ~featSel.(sess).vein_mask & maskROI & featSel.(sess).anyCondActivation_mask;
%             imagesc(maskTmp(:,:,10));
%             sum(maskTmp(:))
        end
        
        % Discriminant voxels
        featSel.(sess).discrim_T2 = nan(size(brain));
        featSel.(sess).discrim_P = nan(size(brain));
        featSel.(sess).discrim_mask = true(size(brain));
        if ~featSel_bSess.vein.doIt
            error('code that')
        else
            % get discriminant voxels (in vector format within ROI)
            x = data.(sess).data;
            ind = false(size(x,[1 2 3]));
            ind(:) = maskROI(maskFitArea);
            x = permute(x,[4 5 1 2 3]);
            x = permute(x(:,:,ind),[1 3 2]); % rep X vox X cond
            y = ones(size(x,1),1);
            x = cat(1,x(:,:,1),x(:,:,2));
            y = cat(1,y.*1,y.*2);
%             spaceList = featSel_bSess.discrim.spaceList;
%             for spaceInd = 1:length(spaceList)
%             switch featSel_bSess.discrim.space
%                 case 'cart'
            nVox = size(x,2);
            x = [real(x) imag(x)];
            featStat = nan(1,nVox);
            featP = nan(1,nVox);
            for voxInd = 1:nVox
                stats = T2Hot2d(x(:,[0 nVox]+voxInd));
                featStat(voxInd) = stats.T2;
                featP(voxInd) = stats.P;
            end
            % convert to map
            featSel.(sess).discrim_T2(maskROI) = featStat;
            featSel.(sess).discrim_P(maskROI) = featP;
            % get thresh based on percentile of ACTIVE NON-VEIN ROI voxels
            mask = ~featSel.(sess).vein_mask & featSel.(sess).anyCondActivation_mask;
            featSel.(sess).discrim_T2thresh = prctile(featStat(mask(maskROI)),20);
            featSel.(sess).discrim_mask = featSel.(sess).discrim_T2>featSel.(sess).discrim_T2thresh;
        end
    end
%     figure('WindowStyle','docked')
%     imagesc(featSel.(sess).anyCondActivation_F(:,:,10))
%     figure('WindowStyle','docked')
%     imagesc(featSel.(sess).anyCondActivation_P(:,:,10))
%     figure('WindowStyle','docked')
%     imagesc(featSel.(sess).anyCondActivation_FDR(:,:,10))
%     figure('WindowStyle','docked')
%     imagesc(featSel.(sess).anyCondActivation_mask(:,:,10))

    
    %% Plot masking
    if subjInd==1
        visibility = 'on';
    elseif figOption.subj==inf && ~verbose
        visibility = 'off';
    else
        visibility = 'on';
    end
    if subjInd==figOption.subj || figOption.subj==inf
        sess = 'sess1';
        slice = 11;
        
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
        im = featSel.(sess).anyCondActivation_F(:,:,slice);
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
        maskTmp = maskROI;
        maskTmp = maskTmp & featSel.(sess).anyCondActivation_mask;
        maskTmp = maskTmp & ~featSel.(sess).vein_mask;
        im(~maskTmp(:,:,slice)) = nan;
        axOver = plotIm(axes,im,cLim);
        alphaData = ~isnan(im);
        makeOverlay(axBak,axOver,alphaData,cMap_F,'log',cLim)
        ylabel({'Activation Level (F value)' ['\fontsize{8}Thresholded (' featSel.(sess).anyCondActivation_threshType '<' num2str(featSel.(sess).anyCondActivation_threshVal,'%0.2f') ') within ROI']});
        % hist
        f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
        vol = featSel.(sess).anyCondActivation_F(maskROI);
        vol = log(vol);
        hHist = histogram(vol,'BinWidth',0.1); hold on
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
        drawnow
        if ~strcmp(featSel.(sess).anyCondActivation_threshType,'none')
            tmp = featSel.(sess).(['anyCondActivation_' (upper(featSel.(sess).anyCondActivation_threshType))])(maskROI);
            [~,b] = min(abs(tmp-featSel.(sess).anyCondActivation_threshVal));
            yLim = ylim;
            plot([1 1].*vol(b),yLim,'r')
            drawnow
            hTex = text(vol(b),yLim(2),['F=' num2str(exp(vol(b)),'%0.1f') ', ' featSel.(sess).anyCondActivation_threshType '=' num2str(featSel.(sess).anyCondActivation_threshVal,'%0.2f')]);
            hTex.VerticalAlignment = 'top';
            hTex.Position(1) = hTex.Position(1) + hTex.Extent(3).*0.02;
            drawnow
        end

        % Vein overlay
        f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
        axBak = plotIm(axes,brain(:,:,slice));
        axOver = plotIm(axes,featSel.(sess).vein_score(:,:,slice));
        alphaData = ~isnan(featSel.(sess).vein_score(:,:,slice));
        cLim = featSel.(sess).vein_score(:,:,slice);
        cLim = cLim(maskROI(:,:,slice));
        cLim = [min(cLim) max(cLim)];
        makeOverlay(axBak,axOver,alphaData,cMap_vein,'lin',cLim)
        ylabel('Veinness (signal std / signal mean)')
        %hist
        f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
        vol = featSel.(sess).vein_score(maskROI);
        hHist = histogram(vol,'BinWidth',0.002); hold on
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
        xlabel('Vein score (BOLD std/mean)')
        drawnow
        if featSel.(sess).vein_doIt
            tmp = featSel.(sess).vein_thresh;
            yLim = ylim;
            plot([1 1].*tmp,yLim,'b')
            hTex = text(tmp,yLim(2),num2str(100-veinPerc,'%0.0f%%ile'));
            hTex.VerticalAlignment = 'top';
            hTex.Position(1) = hTex.Position(1) + hTex.Extent(3).*0.02;
        end

        %V1 retinotopic representation + Vein
        f{end+1} = figure('WindowStyle','docked','color','w','visible',visibility);
        axBak = plotIm(axes,brain(:,:,slice));
        axV1 = plotIm(axes,double(maskROI(:,:,slice)));
        makeOverlay(axBak,axV1,maskROI(:,:,slice),[0 0 0; 255 255 0]./255)
        if featSel.(sess).vein_doIt
            tmp = double(featSel.(sess).vein_mask(:,:,slice));
            tmp(~maskROI(:,:,slice)) = 0;
            axVein = plotIm(axes,tmp);
            makeOverlay(axBak,axVein,tmp,[0 0 0; 0 140 225]./255)
        end

        % Save
        if figOption.save
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
    else
        % clean non-renewed fig
        if figOption.save
            filename = fullfile(pwd,mfilename);
            if ~exist(filename,'dir'); mkdir(filename); end
            filename = fullfile(filename,[subjList{subjInd}]);
            for i = 1:length(f)
                curFile = dir([filename '_' num2str(i) '.*']);
                if ~isempty(curFile)
                    for ii = 1:length(curFile)
                        delFile = fullfile(curFile(ii).folder,curFile(ii).name);
                        if verbose; disp(['delete old: ' delFile]); end
                        delete(delFile)
                    end
                end
            end
        end
    end
    
    
    %% Apply ROI mask and vectorize
    for sessInd = 1:2
        sess = ['sess' num2str(sessInd)];
        sz = size(data.(sess).data);
        voxInd = false(sz(1:3));
        voxInd(:) = maskROI(maskFitArea);
%         figure('WindowStyle','docked')
%         imagesc(voxInd(:,:,10))

        % sin responses
        d.(sess).data = permute(data.(sess).data,[4 5 1 2 3]);
        d.(sess).data = d.(sess).data(:,:,voxInd);
        d.(sess).data = permute(d.(sess).data,[1 3 2]);
        d.(sess).hr = permute(hr.(sess),[4 5 6 1 2 3]);
        d.(sess).hr = d.(sess).hr(:,:,:,voxInd);
        d.(sess).hr = permute(d.(sess).hr,[2 4 3 1]);
        d.(sess).runLabel = permute(data.(sess).runLabel,[4 5 1 2 3]);
        d.(sess).runLabel = permute(d.(sess).runLabel,[1 3 2]);
        d.(sess).info = 'run x vox x cond[ori1, ori2, plaid] X TR';

        % feature slection stats
        fieldList = fields(featSel.(sess));
        for i = 1:length(fieldList)
            isDataField = (isnumeric(featSel.(sess).(fieldList{i})) || islogical(featSel.(sess).(fieldList{i})) ) && length(size(featSel.(sess).(fieldList{i})))==3 && all(size(featSel.(sess).(fieldList{i}))==size(brain));
            if isDataField
                tmp = nan(sz(1:3));
                tmp(:) = featSel.(sess).(fieldList{i})(maskFitArea);
                featSel.(sess).(fieldList{i}) = tmp(voxInd)';
            end
        end
        featSel.(sess).info = '1 x vox';
    end
    clear hr

    %% Simplify stuff
    for sessInd = 1:2
        sess = ['sess' num2str(sessInd)];
        fieldList = fields(featSel.(sess));
        for i = 1:length(fieldList)
            d.(sess).(fieldList{i}) = featSel.(sess).(fieldList{i});
        end
    end
    clear featSel

    %% Export some parameters
    param.subjList = subjList;
    param.brain = brain;
    param.anyCondActivation.doIt = d.(sess).anyCondActivation_doIt;
    param.anyCondActivation.fitType = d.(sess).anyCondActivation_fitType;
    param.anyCondActivation.threshType = d.(sess).anyCondActivation_threshType;
    param.anyCondActivation.threshVal = d.(sess).anyCondActivation_threshVal;
    param.vein.doIt = d.(sess).vein_doIt;
    param.vein.source = d.(sess).vein_source;
    param.vein.perc = d.(sess).vein_perc;

    %% Save
    if ~exist(fullfile(funPath,funLevel3),'dir')
        mkdir(fullfile(funPath,funLevel3));
    end

    tmp = fullfile(funPath,funLevel3,[subjList{subjInd} '_' mfilename]);
    if verbose; disp(['Saving to: ' tmp '.mat']); end
    save(tmp,'d','param')
end
