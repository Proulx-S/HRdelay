function maskSinAndHrFit(fitType,threshType)
noMovement = 1;
actuallyRun = 1;
saveFig = 1;
plotAllSubj = 1;
doVein = 1;
veinSource = 'reducedModelResid'; % 'reducedModelResid' (stimulus-driven signal included in std) or 'fullModelResid (stimulus-driven signal excluded in std)'
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
    if doVein
        sessLabel = [results.inputs.opt.sessionLabel{:}]';
        for sessInd = 1:2
            vein.(['sess' num2str(sessInd)]).noiseOverMean = nan(size(maskAnat));
            switch veinSource
                case 'reducedModelResid'
                    vein.(['sess' num2str(sessInd)]).noiseOverMean(maskFit) = mean(results.OLS.mixed.veinFull(:,:,:,sessLabel==sessInd),4);
                case 'fullModelResid'
                    vein.(['sess' num2str(sessInd)]).noiseOverMean(maskFit) = mean(results.OLS.mixed.veinReduced(:,:,:,sessLabel==sessInd),4);
            end
            vein.(['sess' num2str(sessInd)]).thresh = prctile(vein.(['sess' num2str(sessInd)]).noiseOverMean(maskAnat),100-veinPerc);
            vein.(['sess' num2str(sessInd)]).mask = vein.(['sess' num2str(sessInd)]).noiseOverMean>vein.(['sess' num2str(sessInd)]).thresh;
        end
    end
    
    % feature selection
    for sessInd = 1:2
        switch fitType
            case 'mixed'
                F = results.OLS.(fitType).(['Fsess' num2str(sessInd)]).val.F;
                df = results.OLS.(fitType).(['Fsess' num2str(sessInd)]).df;
            case 'fixed'
                F = results.OLS.(fitType).(['sess' num2str(sessInd)]).F.val.F;
                df = results.OLS.(fitType).(['sess' num2str(sessInd)]).F.df;
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
%         cMap_F = brewermap(256,'reds');
%         cMap_vein = brewermap(256,'blues');
%         save(filename,'cMap_F','cMap_vein');
        load(filename,'cMap_F','cMap_vein');
        
        f = [];
        %Brain
        f(end+1) = figure('WindowStyle','docked','color','w');
        axBak = plotIm(axes,brain(:,:,slice));
        title('BOLD image (1 TR)')
        
        %Activation overlay
        f(end+1) = figure('WindowStyle','docked','color','w');
        axBak = plotIm(axes,brain(:,:,slice));
        im = featSel.(sess).anyCondActivation.F(:,:,slice);
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
        maskThresh = featSel.(sess).anyCondActivation.(upper(threshType))<threshVal;
        if doVein
            maskVein = vein.(sess).mask;
            maskFull = maskThresh & ~maskVein & maskAnat;
        else
            maskFull = maskThresh & maskAnat;
        end
        im(~maskFull(:,:,slice)) = nan;
        axOver = plotIm(axes,im,cLim);
        alphaData = ~isnan(im);
        makeOverlay(axBak,axOver,alphaData,cMap_F,'log',cLim)
        ylabel({'Activation Level (F value)' ['\fontsize{8}Thresholded (' threshType '<' num2str(threshVal,'%0.2f') ') within ROI']});
        %hist
        f(end+1) = figure('WindowStyle','docked','color','w');
        vol = featSel.(sess).anyCondActivation.F(maskAnat);
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
        tmp = featSel.(sess).anyCondActivation.(upper(threshType))(maskAnat);
        [~,b] = min(abs(tmp-threshVal));
        yLim = ylim;
        plot([1 1].*vol(b),yLim,'r')
        hTex = text(vol(b),yLim(2),['F=' num2str(exp(vol(b)),'%0.1f') ', ' threshType '=' num2str(threshVal,'%0.2f')]);
        hTex.VerticalAlignment = 'top';
        hTex.Position(1) = hTex.Position(1) + hTex.Extent(3).*0.02;
        
        
        if doVein
            %Vein overlay
            f(end+1) = figure('WindowStyle','docked','color','w');
            axBak = plotIm(axes,brain(:,:,slice));
            axOver = plotIm(axes,vein.(sess).noiseOverMean(:,:,slice));
            alphaData = ~isnan(featSel.(sess).anyCondActivation.F(:,:,slice));
            cLim = vein.(sess).noiseOverMean(:,:,slice);
            cLim = cLim(maskAnat(:,:,slice));
            cLim = [min(cLim) max(cLim)];
            makeOverlay(axBak,axOver,alphaData,cMap_vein,'lin',cLim)
            ylabel('Veinness (signal std / signal mean)')
            %hist
            f(end+1) = figure('WindowStyle','docked','color','w');
            vol = vein.(sess).noiseOverMean(maskAnat);
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
        axV1 = plotIm(axes,double(maskAnat(:,:,slice)));
        makeOverlay(axBak,axV1,maskAnat(:,:,slice),[0 0 0; 255 255 0]./255)
        tmp = double(vein.(sess).mask(:,:,slice));
        tmp(~maskAnat(:,:,slice)) = 0;
        axVein = plotIm(axes,tmp);
        makeOverlay(axBak,axVein,tmp,[0 0 0; 0 140 225]./255)
        
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


