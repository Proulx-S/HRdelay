function maskSinFit
actuallyRun = 1;
plotAll = 0;
threshType = 'p'; % 'none', 'p' or 'fdr'
threshVal = 0.05;

%% Define paths
repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
dataDir = 'C-derived\DecodingHR';
anatPath = fullfile(repoPath,dataDir,'anat');
anatLevel = 'z';
funPath = fullfile(repoPath,dataDir,'fun');
funLevel1 = 'x';
funLevel2 = 'y';
funLevel3 = 'zSin';
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';

if ~actuallyRun
    disp('Not actually running to save some time')
end
disp(['IN: anatomical V1 roi (C-derived\DecodingHR\anat\' anatLevel ')'])
disp(['IN: voxel visual field eccentricity (C-derived\DecodingHR\anat\' anatLevel ')'])
disp(['IN: sinusoidal fit results (C-derived\DecodingHR\fun\' funLevel2 ')'])
disp('F(IN)=OUT: masks the fit according to voxel eccentricity and activation level')
disp(['OUT: sinusoidal fit results (C-derived\DecodingHR\fun\' funLevel3 ')'])


for subjInd = 1:length(subjList)
    %% Load fun data
    load(fullfile(funPath,funLevel2,subjList{subjInd},'v1SinCos_1perRun_move12.mat'),'results')
    % brain
    tmpFilename = ls(fullfile(funPath,funLevel1,subjList{subjInd},'trun101_preprocessed.nii*'));
    a = load_nii(fullfile(funPath,funLevel1,subjList{subjInd},tmpFilename));
    brain = flip(permute(a.img,[3 1 2 4]),1); clear a
    brain = brain(:,:,:,1);
    
    %% Set masks
    % maskV1
    tmpFilename = ls(fullfile(anatPath,anatLevel,subjList{subjInd},'v1.nii*'));
    a = load_nii(fullfile(anatPath,anatLevel,subjList{subjInd},tmpFilename));
    maskV1 = flipdim(permute(a.img,[3 1 2 4]),1); clear a
    maskV1(:,:,1) = zeros(size(maskV1,1),size(maskV1,2));% Remove corrupted slices
    maskV1(:,:,end) = zeros(size(maskV1,1),size(maskV1,2));% Remove corrupted slices
    
    % maskECC
    tmpFilename = ls(fullfile(anatPath,anatLevel,subjList{subjInd},'lh.ecc.nii*'));
    a = load_nii(fullfile(anatPath,anatLevel,subjList{subjInd},tmpFilename));
    eccL = flipdim(permute(a.img,[3 1 2 4]),1); clear a
    tmpFilename = ls(fullfile(anatPath,anatLevel,subjList{subjInd},'rh.ecc.nii*'));
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
        eval(['F = results.OLS.mixed.Fsess' num2str(sessInd) '.val.F(maskAnat(maskFit));']);
        eval(['df = results.OLS.mixed.Fsess' num2str(sessInd) '.df;']);
        
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
    if plotAll || subjInd==1
        sess = 'sess1';
        
        f = figure('WindowStyle','docked');
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
        switch threshType
            case 'p'
                thresh = maskFun.(sess).P(:,:,10);
            case 'fdr'
                thresh = maskFun.(sess).FDR(:,:,10);
            otherwise
                error('x')
        end
        F(thresh>threshVal) = 0; % functional threshold here
        F(isnan(F)) = 0; % functional threshold here
        im2 = ind2rgb(uint8(F./max(F(:))*255),autumn(256));
        ind = repmat(F~=0,[1 1 3]);
        im1(ind) = im2(ind);
        imshow(im1)
        ax3.PlotBoxAspectRatio = [1 1 1];
        
        axCB = copyobj(ax3,gcf);
        axes(axCB);
        CB = permute(flipud(autumn(size(ax1.Children.CData,1))),[1 3 2]);
        CB = repmat(CB,[1 round(0.1*size(ax1.Children.CData,2)) 1]);
        imshow(CB);
        axCB.Position(1) = ax3.Position(1) - ax3.Position(3)*0.5;
        axCB.YAxis.Visible = 'on';
        axCB.YAxis.TickValues = [1 size(ax1.Children.CData,1)];
        tmp = threshVal-thresh(:); tmp(tmp<0) = nan; [~,b] = min(tmp);
        axCB.YAxis.TickLabels = {num2str(max(F(:)),'%0.1f') num2str(F(b),'%0.1f')};
        axCB.YAxis.Label.String = {'F' ['p<' num2str(threshVal,'%0.2f')]};
        axCB.YAxis.Label.Rotation = 0;
        axCB.YAxis.Label.VerticalAlignment = 'middle';
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
        line([1 1].*F(b),ylim,'Color','r');
        text(F(b)+diff(xlim)*0.01,ax4.YLim(2)*0.95,['p=' num2str(threshVal,'%0.2f')])
        title('All ROI voxels')
        ax4.XTick = sort([ax4.XTick round(F(b),2)]);
        ax4.TickDir = 'out';
        ax4.XTick(ax4.XTick==0) = [];
        
        suptitle([subjList{subjInd} '; ' sess])
        f.Color = 'none';
        drawnow
    end
    
    %% Stop here if not actually runing
    if ~actuallyRun
        return
    end
    
    %% Apply masks to fun data, vectorize voxels and split sessions and conditions
    [X,Y] = pol2cart(results.OLS.mixed.delay,results.OLS.mixed.amp);
    data = complex(X,Y); clear X Y
    
    % mask and vectorize
    sz = size(data);
    dataX = nan(sz(4),sum(maskAnat(maskFit)));
    for runInd = 1:sz(4)
        curData = data(:,:,:,runInd);
        tmp = false(size(curData));
        tmp(:) = maskAnat(maskFit);
        dataX(runInd,:) = curData(tmp);
    end
    list = {'F' 'FDR' 'P'};
    for i = 1:3
        for sessInd = 1:2
            tmp = maskFun.(['sess' num2str(sessInd)]).(list{i});
            tmp = tmp(maskFit);
            tmp = tmp(maskAnat(maskFit))';
            stats.(list{i}).(['sess' num2str(sessInd)]) = tmp;
        end
    end
    
    % split sessions and conditions
    sessLabel = [results.inputs.opt.sessionLabel{:}];
    condLabel = repmat(1:3,[size(dataX,1)/3 1]); condLabel = condLabel(:)';
    for sessInd = 1:2
        d.(['sess' num2str(sessInd)]).xData = nan(sum(sessLabel==sessInd)/3,size(dataX,2),3);
        for condInd = 1:3
            d.(['sess' num2str(sessInd)]).xData(:,:,condInd) = dataX(sessLabel==sessInd & condLabel==condInd,:);
        end
        d.(['sess' num2str(sessInd)]).info = '%BOLD: run x vox x cond[ori1, ori2, plaid]';
        for i = 1:3
            d.(['sess' num2str(sessInd)]).(list{i}) = stats.(list{i}).(['sess' num2str(sessInd)]);
        end
    end
    
    
    % new = cat(1,d.(sess).xData(:,:,1),d.(sess).xData(:,:,2),d.(sess).xData(:,:,3));
    % old = load(['C:\Users\sebas\OneDrive - McGill University\McGill\work\projects\170210_HRdecoding\C_processing\02jp_' sess '.mat']);
    % old = cat(1,old.d.xData,old.d.normData);
    % scatter(abs(old(:)),abs(new(:))); hold on
    % line(xlim,ylim)
    
    
    %% Save
    if ~exist(fullfile(funPath,funLevel3),'dir')
        mkdir(fullfile(funPath,funLevel3));
    end
    
    save(fullfile(funPath,funLevel3,[subjList{subjInd} '_' mfilename]),'d')
end


