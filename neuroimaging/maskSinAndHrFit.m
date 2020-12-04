function maskSinAndHrFit(fitType,threshType)
noMovement = 1;
actuallyRun = 1;
saveFig = 0;
plotAllSubj = 0;
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
    
    % maskFun
    for sessInd = 1:2
        switch fitType
            case 'mixed'
                eval(['F = results.OLS.mixed.Fsess' num2str(sessInd) '.val.F(maskAnat(maskFit));']);
                eval(['df = results.OLS.mixed.Fsess' num2str(sessInd) '.df;']);
            case 'fixed'
                eval(['F = results.OLS.fixed.sess' num2str(sessInd) '.F.val.F(maskAnat(maskFit));']);
                eval(['df = results.OLS.fixed.sess' num2str(sessInd) '.F.df;']);
            otherwise
                error('X')
        end
        
        [FDR,P] = getPfromF(F,df);
        tmp = nan(size(brain)); tmp(maskAnat) = F;
        F = tmp;
        
        tmp = nan(size(brain)); tmp(maskAnat) = FDR;
        FDR = tmp;
        
        tmp = nan(size(brain)); tmp(maskAnat) = P;
        P = tmp;
        
        maskFun.(['sess' num2str(sessInd)]).F = F;
        maskFun.(['sess' num2str(sessInd)]).P = P;
        maskFun.(['sess' num2str(sessInd)]).FDR = FDR;
    end
    
    %% Plot masking
    if plotAllSubj || subjInd==1
        sess = 'sess1';
        
        fSubj(subjInd) = figure('WindowStyle','docked');
        ax1 = subplot(2,2,1);
        imagesc(brain(:,:,10)); colormap gray; axis off;
        ax1.PlotBoxAspectRatio = [1 1 1];
        title('BOLD image (1 TR)')
        
        ax2 = subplot(2,2,2);
        imagesc(maskAnat(:,:,10)); colormap gray; axis off
        ax2.PlotBoxAspectRatio = [1 1 1];
        title({'Anatomical ROI' '\fontsize{8} V1, 1 to 6 dva eccentricity'})
        
        ax3 = subplot(2,2,3); axis off
        im1 = brain(:,:,10);
        im1 = ind2rgb(uint8(im1./max(im1(:))*255),gray(256));
        F = maskFun.(sess).F(:,:,10);
        F(isnan(F)) = 0;
        switch threshType
            case 'p'
                thresh = maskFun.(sess).P(:,:,10);
            case 'fdr'
                thresh = maskFun.(sess).FDR(:,:,10);
            case 'none'
            otherwise
                error('x')
        end
        if ~strcmp(threshType,'none')
            F(thresh>threshVal) = 0; % functional threshold here
        end
        tmp = unique(F(:));
        [~,b] = sort(tmp);
        Fmin =  tmp(b(2)); clear b
        Fmax = max(F(:));
        im2 = ind2rgb(uint8((F-Fmin)./(Fmax-Fmin)*255),autumn(256));
        ind = repmat(F~=0,[1 1 3]);
        im1(ind) = im2(ind);
        imshow(im1);
        ax3.PlotBoxAspectRatio = [1 1 1];
        
        axCB = copyobj(ax3,gcf);
        axes(axCB);
        CB = permute(flipud(autumn(size(ax1.Children.CData,1))),[1 3 2]);
        CB = repmat(CB,[1 round(0.1*size(ax1.Children.CData,2)) 1]);
        imshow(CB);
        axCB.Position(1) = ax3.Position(1) - ax3.Position(3)*0.5;
        axCB.YAxis.Visible = 'on';
        axCB.YAxis.TickValues = [1 size(ax1.Children.CData,1)];
        if strcmp(threshType,'none')
            axCB.YAxis.Label.String = 'F';
            axCB.YAxis.TickLabels = {num2str(Fmax,'%0.1f') num2str(Fmin,'%0.1f')};
        else
            tmp = threshVal-thresh(:); tmp(tmp<0) = nan; [~,b] = min(tmp);
            axCB.YAxis.TickLabels = {num2str(max(F(:)),'%0.1f') num2str(F(b),'%0.1f')};
            if strcmp(threshType,'p')
                axCB.YAxis.Label.String = {'F' ['p<' num2str(threshVal,'%0.2f')]};
            elseif strcmp(threshType,'fdr')
                axCB.YAxis.Label.String = {'F' ['FDR<' num2str(threshVal,'%0.2f')]};
            end
        end
        axCB.YAxis.Label.Rotation = 90;
        axCB.YAxis.Label.VerticalAlignment = 'top';
        axCB.YAxis.Label.HorizontalAlignment = 'center';
        
        ax3.Title.String = 'Visual Activation';
        
        ax4 = subplot(2,2,4);
        hm = histogram(maskFun.(sess).F(maskAnat),1000); hold on
        xlim([0 20])
        ylim(ylim*1.05)
        hm.FaceColor = 'k';
        ax4.PlotBoxAspectRatio = [1 1 1];
        ylabel('vox count');
        xlabel('F')
        switch threshType
            case 'none'
            case 'p'
                line([1 1].*F(b),ylim,'Color','r');
                text(F(b)+diff(xlim)*0.01,ax4.YLim(2)*0.95,['p=' num2str(threshVal,'%0.2f')])
            case 'fdr'
                line([1 1].*F(b),ylim,'Color','r');
                text(F(b)+diff(xlim)*0.01,ax4.YLim(2)*0.95,['FDR=' num2str(threshVal,'%0.2f')])
            otherwise
                error('X')
        end
        title('All ROI voxels')
        if ~strcmp(threshType,'none')
            ax4.XTick = sort([ax4.XTick round(F(b),2)]);
        end
        ax4.TickDir = 'out';
        ax4.XTick(ax4.XTick==0) = [];
        
        suptitle([subjList{subjInd} '; ' sess])
        
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
    [X,Y] = pol2cart(results.OLS.mixed.delay,results.OLS.mixed.amp);
    data = complex(X,Y); clear X Y
    
    hr.sess1 = resultsResp.OLS.mixed.sess1.resp;
    hr.sess2 = resultsResp.OLS.mixed.sess2.resp;
    hr.info = 'x X y X z X TR X run X cond[ori1,ori2,plaid]';
    
    % mask and vectorize
    % sin responses
    sz = size(data);
    dataX = nan(sz(4),sum(maskAnat(maskFit)));
    runLabelX = nan(sz(4),1);
    for runInd = 1:sz(4)
        curData = data(:,:,:,runInd);
        tmp = false(size(curData));
        tmp(:) = maskAnat(maskFit);
        dataX(runInd,:) = curData(tmp);
        runLabelX(runInd,1) = results.inputs.opt.runLabel{runInd};
    end
    % stats
    list = {'F' 'FDR' 'P'};
    for i = 1:3
        for sessInd = 1:2
            tmp = maskFun.(['sess' num2str(sessInd)]).(list{i});
            tmp = tmp(maskFit);
            tmp = tmp(maskAnat(maskFit))';
            stats.(list{i}).(['sess' num2str(sessInd)]) = tmp;
        end
    end
    % response shape
    for sessInd = 1:2
        tmp = permute(hr.(['sess' num2str(sessInd)]),[4 5 6 1 2 3]);
        tmpMask = false(size(tmp,4:6));
        tmpMask(:) = maskAnat(maskFit);
        hr.(['sess' num2str(sessInd)]) = permute(tmp(:,:,:,tmpMask),[2 4 3 1]);
        clear tmp tmpMask
    end
    hr.info = '%BOLD: run X vox X cond[ori1,ori2,plaid] X TR';
%   
    
    % split sessions and conditions
    sessLabel = [results.inputs.opt.sessionLabel{:}];
    condLabel = repmat(1:3,[size(dataX,1)/3 1]); condLabel = condLabel(:)';
    for sessInd = 1:2
        d.(['sess' num2str(sessInd)]).xData = nan(sum(sessLabel==sessInd)/3,size(dataX,2),3);
        for condInd = 1:3
            d.(['sess' num2str(sessInd)]).xData(:,:,condInd) = dataX(sessLabel==sessInd & condLabel==condInd,:);
        end
        d.(['sess' num2str(sessInd)]).info = '%BOLD: run x vox x cond[ori1, ori2, plaid]';
        d.(['sess' num2str(sessInd)]).runLabel = nan(sum(sessLabel==sessInd)/3,1,3);
        for condInd = 1:3
            d.(['sess' num2str(sessInd)]).runLabel(:,:,condInd) = runLabelX(sessLabel==sessInd & condLabel==condInd,1);
        end
        for i = 1:3
            d.(['sess' num2str(sessInd)]).(list{i}) = stats.(list{i}).(['sess' num2str(sessInd)]);
        end
    end
    d.sess2.runLabel = d.sess2.runLabel - min(d.sess2.runLabel(:)) + 1;
    
    
    % new = cat(1,d.(sess).xData(:,:,1),d.(sess).xData(:,:,2),d.(sess).xData(:,:,3));
    % old = load(['C:\Users\sebas\OneDrive - McGill University\McGill\work\projects\170210_HRdecoding\C_processing\02jp_' sess '.mat']);
    % old = cat(1,old.d.xData,old.d.normData);
    % scatter(abs(old(:)),abs(new(:))); hold on
    % line(xlim,ylim)
    
    
    %% Save
    if ~exist(fullfile(funPath,funLevel3Sin),'dir')
        mkdir(fullfile(funPath,funLevel3Sin));
    end
    if ~exist(fullfile(funPath,funLevel3Hr),'dir')
        mkdir(fullfile(funPath,funLevel3Hr));
    end
    
    if noMovement
        save(fullfile(funPath,funLevel3Sin,[subjList{subjInd} '_' mfilename '_noMovement']),'d')
        save(fullfile(funPath,funLevel3Hr,[subjList{subjInd} '_' mfilename '_noMovement']),'hr')
    else
        save(fullfile(funPath,funLevel3Sin,[subjList{subjInd} '_' mfilename]),'d')
        save(fullfile(funPath,funLevel3Hr,[subjList{subjInd} '_' mfilename]),'hr')
    end
end


