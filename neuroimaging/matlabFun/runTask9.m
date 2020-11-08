function [ taskDone ] = runTask9(allFile,perm,altMainDir,guillimin,ppn)
global tmpStrResp
tmpStrResp = {'respSVM__pcaNone_pcaNone','respSVM__pcaRun_pcaNone','respSVM__pcaTime_pcaNone','respSVM__pcaNone_pcaSpace','respSVM__pcaRun_pcaSpace','respSVM__pcaTime_pcaSpace','respSVM__pcaNone_pcaRun','respSVM__pcaRun_pcaRun','respSVM__pcaTime_pcaRun','respPCA','respPCA_PCAalign'};
           
addpath(genpath('C:\Users\Sebastien\OneDrive - McGill University\MATLAB\libsvm-3.21'))

rng('shuffle')
if ~exist('guillimin','var')
    guillimin=0;
end
if ~exist('ppn','var')
    ppn=0;
end

taskDone = false(size(allFile));
for taskInd = 1:length(allFile)
    clearvars -except taskInd allFile plotIt doPerm dataDirOut perm altMainDir guillimin ppn taskDone tmpStrResp
    display(['Performing ' char(allFile{taskInd})])
    
    if length(find(size(allFile)>1))>1
        load(allFile{taskInd,1})
        load(allFile{taskInd,2})
    else
        load(allFile{taskInd})
    end
    
    
    if guillimin
        dataDirIn = fullfile(remoteDir,subj);
        dataDirOut = remoteDir;
    end
    
    %% Do SVM if not already done
    stopIt = 0;
    stopThis = 0;
    if ~exist('svm','var')
        stopIt = 1;
        
        %% Load data, compile it and threshold it
        %Compile it
        switch fitType
            case 'SinCos'
                %Load mask (and brain for reference)
                a = load_nii(fullfile(dataDirIn,'run1c_masking',ls(fullfile(dataDirIn,'run1c_masking',[p.maskName '.nii*']))));
                metric.mask = int8(flipdim(permute(a.img,[3 1 2 4]),1)); clear a
                metric.mask(:,:,1) = zeros(size(metric.mask,1),size(metric.mask,2));% Remove corrupted slices
                metric.mask(:,:,end) = zeros(size(metric.mask,1),size(metric.mask,2));% Remove corrupted slices
                a = load_nii(fullfile(dataDirIn,'run1a_preprocessing',ls(fullfile(dataDirIn,'run1a_preprocessing','trun101_preprocessed.nii*'))));
                metric.brain = mean(flipdim(permute(a.img,[3 1 2 4]),1),4); clear a
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Should review that if want to include ecc
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Apply F thresh to mask
%                 load(fullfile(dataDirIn,'sepSess170108_v1SinCos_1perRun_move12.mat'),'results')
%                 load(fullfile(dataDirIn,'sepSess170108_v1SinCos_1perRun_move12.mat'),'results')
                load(fullfile(dataDirIn,['170505_v1' fitType '_' p.split '_' p.motionParam p.sm '.mat']),'results')
                fitMask = results.mask;
                anatMask = logical(metric.mask);
                anatMask = anatMask(:,:,squeeze(any(any(fitMask,1),2)),:);
                anatMask = anatMask(:,squeeze(any(any(fitMask,1),3)),:,:);
                anatMask = anatMask(squeeze(any(any(fitMask,2),3)),:,:,:);
                
                Fana = 'fixedAll3';
%                 Fana = 'mixedAll3';
%                 Fana = 'mixedCond3';
                switch Fana
                    case 'fixedAll3'
                        F = results.OLS.fixed.(p.preSelection).F.val.F;
                        [FDR,~] = getPfromF(F(logical(anatMask)),results.OLS.fixed.(p.preSelection).F.df);
                    case 'mixedAll3'
                        F = results.OLS.mixed.(['F' p.preSelection]).val.F;
                        [FDR,~] = getPfromF(F(logical(anatMask)),results.OLS.mixed.(['F' p.preSelection]).df);
                    case 'mixedCond3'
                        F = results.OLS.mixed.(['Fcond3' p.preSelection]).val.F;
                        [FDR,~] = getPfromF(F(logical(anatMask)),results.OLS.mixed.(['Fcond3' p.preSelection]).df);
                end
                
                %convert from vector to fit volume
                tmpFDR1 = nan(size(F));
                tmpFDR1(logical(anatMask)) = FDR; FDR = tmpFDR1;
                %convert from fit volume to imaging volume
                tmpFDR2 = nan(size(fitMask));
                tmpFDR2(fitMask) = FDR; FDR = tmpFDR2;
                
                
                Fmask = FDR>0.05;
                if p.thresholdData
                    metric.mask(Fmask) = 0;
                end
                
                
                fixedF = results.OLS.fixed.(p.preSelection).F.val.F;
                [FDR,~] = getPfromF(F(logical(anatMask)),results.OLS.fixed.(p.preSelection).F.df);
                fixedFDR = nan(size(F));
                fixedFDR(logical(anatMask)) = FDR;
                %%% Get fixed effect fit
                fixedAmp = results.OLS.fixed.(p.preSelection).amp;
                fixedDelay = results.OLS.fixed.(p.preSelection).delay;
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                
                %Load ecc mask
                a = load_nii(fullfile(dataDirIn,'run1c_masking',['lh.ecc.nii.gz']));
                tmpL = flipdim(permute(a.img,[3 1 2 4]),1); clear a
                a = load_nii(fullfile(dataDirIn,'run1c_masking',['rh.ecc.nii.gz']));
                tmpR = flipdim(permute(a.img,[3 1 2 4]),1); clear a
                tmpL(find(tmpL~=0 & tmpR~=0))=0;
                metric.ecc = tmpL+tmpR; clear tmpL tmpR
                metric.ecc(:,:,1) = zeros(size(metric.mask,1),size(metric.mask,2));% Remove corrupted slices
                metric.ecc(:,:,end) = zeros(size(metric.mask,1),size(metric.mask,2));% Remove corrupted slices
                
                %Load pol mask
                a = load_nii(fullfile(dataDirIn,'run1c_masking',['lh.pol.nii.gz']));
                tmpL = flipdim(permute(a.img,[3 1 2 4]),1); clear a
                a = load_nii(fullfile(dataDirIn,'run1c_masking',['rh.pol.nii.gz']));
                tmpR = flipdim(permute(a.img,[3 1 2 4]),1); clear a
                tmpR(tmpR==0)=nan;
                tmpR = (180-tmpR)+180;
                tmpR(isnan(tmpR)) = 0;
                tmpL(find(tmpL~=0 & tmpR~=0))=0;
                metric.pol = tmpL+tmpR; clear tmpL tmpR
                metric.pol(:,:,1) = zeros(size(metric.mask,1),size(metric.mask,2));% Remove corrupted slices
                metric.pol(:,:,end) = zeros(size(metric.mask,1),size(metric.mask,2));% Remove corrupted slices
                
                if p.useEccMask
                    eccMask = false(size(metric.mask));
                    eccMask(metric.ecc>=p.eccLowLim&metric.ecc<p.eccHighLim) = true;
                    metric.mask = metric.mask.*int8(eccMask); clear eccMask
                end
                
                %Load data of interest
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 load(fullfile(dataDirIn,['161212_' upper(p.maskName) fitType '_' p.split '_' p.motionParam p.sm '.mat']),'results')
                %                 load(fullfile(dataDirIn,['161212_v1v2v3' fitType '_' p.split '_' p.motionParam p.sm '.mat']),'results')
                %                 load(fullfile(dataDirIn,['X161212_v1' fitType '_' p.split '_' p.motionParam p.sm '.mat']),'results')
                load(fullfile(dataDirIn,['170505_v1' fitType '_' p.split '_' p.motionParam p.sm '.mat']),'results')
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                respFlag = 0;
                respFlagFix = 0;
                respFlagMix = 0;
                
                for ii = 1:length(p.infoComb)
                    switch p.infoComb{ii}
                        case {'respSVM__pcaNone_pcaNone' 'respSVM__pcaRun_pcaNone' 'respSVM__pcaTime_pcaNone' 'respSVM__pcaNone_pcaSpace' 'respSVM__pcaRun_pcaSpace' 'respSVM__pcaTime_pcaSpace' 'respSVM__pcaNone_pcaRun' 'respSVM__pcaRun_pcaRun' 'respSVM__pcaTime_pcaRun'}
                            respFlag = 1;
                        case {'pcaOnResp_lowLevel' 'pcaOnResp_pcaOnSpace'}
                            respFlag = 1;
                        case {'respPCA','respPCA_PCAalign'}
                            respFlag = 1;
                        case {'pcaOnTS_lowLevel','lowLevel','crossValDelay','crossValDelayQ'}
                        otherwise
                            error('XX');
                    end
                end
                if isfield(p,'saveDataOnly') && p.saveDataOnly
                    respFlagFix = 1;
                    respFlagMix = 1;
                    respFlag = 1;
                end
                respFlag2 = 0;
                for ii = 1:length(p.infoComb)
                    switch p.infoComb{ii}
                        case 'pcaOnTS_lowLevel'
                            respFlag2 = 1;
                        case {'respSVM__pcaNone_pcaNone' 'respSVM__pcaRun_pcaNone' 'respSVM__pcaTime_pcaNone' 'respSVM__pcaNone_pcaSpace' 'respSVM__pcaRun_pcaSpace' 'respSVM__pcaTime_pcaSpace' 'respSVM__pcaNone_pcaRun' 'respSVM__pcaRun_pcaRun' 'respSVM__pcaTime_pcaRun'}
                        case {'pcaOnResp_lowLevel' 'pcaOnResp_pcaOnSpace','lowLevel'}
                        case {'respPCA','respPCA_PCAalign','crossValDelay','crossValDelayQ'}
                        otherwise
                            error('XX');
                    end
                end
                if isfield(p,'saveDataOnly') && p.saveDataOnly
                    respFlag2 = 0;
                end
                if respFlagFix || respFlagMix
                    tmp = results;
                    load(fullfile(dataDirIn,['170505_v1resp_' p.split '_' p.motionParam p.sm '_resp.mat']),'results')
                    if respFlagFix
                        respFix1 = results.OLS.fixed.sess1.resp;
                        respFix2 = results.OLS.fixed.sess2.resp;
                    end
                    if respFlagMix
                        respMix1 = results.OLS.mixed.sess1.resp;
                        respMix2 = results.OLS.mixed.sess2.resp;
                    end
                    clear results; results = tmp; clear tmp
                end
                if respFlag || respFlag2
                    load(fullfile(dataDirIn,['170505_v1' fitType '_' p.split '_' p.motionParam p.sm '_detrend.mat']),'dataDetrend')
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if find(fitMask~=results.mask)
                    error('something wrong with the fitMask used for Fthresh and the the fitMask used for SVM data')
                end
                
                %Define some flags to extract data correctly
                ind = find(fitMask);
                switch p.preSelection
                    case 'sess1'
                        Ffield2 = 'Fcond3sess1';
                        if respFlagFix
                            respFix = respFix1; clear respFix1 respFix2
                        end
                        if respFlagMix
                            respMix = respMix1; clear respMix1 respMix2
                        end
                    case 'sess2'
                        Ffield2 = 'Fcond3sess2';
                        if respFlagFix
                            respFix = respFix2; clear respFix1 respFix2
                        end
                        if respFlagMix
                            respMix = respMix2; clear respMix1 respMix2
                        end
                    case 'bothSess'
                        Ffield2 = [];
                    otherwise
                        error('did not program this yet')
                end
                nRuns = length(results.inputs.datasize)/3;
                nPtsPerRun = results.inputs.datasize{1}(4);
                
                %Get F values (and put it in full VOI)
                if p.regSession
                    Ffield1 = 'mixed_sessReg';
                else
                    Ffield1 = 'mixed';
                end
                
                metric.allExp.F = nan(size(metric.mask));
                %                 if ~isempty(Ffield2)
                %                     metric.allExp.F(ind) = results.OLS.(Ffield1).(Ffield2).val.F;
                %                 end
                
                
                %Get amp and delay (and put it in full VOI)
                if p.regSession
                    Ffield1 = 'mixed_sessReg';
                else
                    Ffield1 = 'mixed';
                end
                
                if respFlagFix
                    metric.respFix = nan([size(metric.mask) 12 3]); % X x Y x Z x t x cond x run
                end
                if respFlagMix
                    metric.respMix = nan([size(metric.mask) 12 3 size(respMix,5)]); % X x Y x Z x t x cond x run
                end
                if respFlag
                    metric.resp = nan([size(metric.mask) 12 3 nRuns]); % X x Y x Z x t x cond x run
                end
                if respFlag2
                    metric.resp2 = nan([size(metric.mask) 12*8 3 nRuns]); % X x Y x Z x t x cond x run
                end
%                 if respFlagFix || respFlagMix
%                     sz = size(respFix);
%                     sz = sz(1:3);
                if respFlag || respFlag2
                    sz = size(dataDetrend{1});
                end
                for cond = 1:3
                    %mixed-effect
                    tmpAmp = results.OLS.(Ffield1).amp(:,:,:,1+(cond-1)*nRuns:cond*nRuns);
                    tmpDelay = results.OLS.(Ffield1).delay(:,:,:,1+(cond-1)*nRuns:cond*nRuns);
                    for run = 1:nRuns
                        tmp = nan(size(metric.mask));
                        tmp(ind) = tmpAmp(:,:,:,run);
                        metric.allRuns.amp(:,:,:,run,cond) = tmp;
                        tmp = nan(size(metric.mask));
                        tmp(ind) = tmpDelay(:,:,:,run);
                        metric.allRuns.delay(:,:,:,run,cond) = tmp;
                    end
                    %fixed-effect
                    tmpAmp = fixedAmp(:,:,:,cond);
                    tmpDelay = fixedDelay(:,:,:,cond);
                    tmp = nan(size(metric.mask));
                    tmp(ind) = tmpAmp;
                    metric.allExp.amp(:,:,:,1,cond) = tmp;
                    tmp = nan(size(metric.mask));
                    tmp(ind) = tmpDelay;
                    metric.allExp.delay(:,:,:,1,cond) = tmp;
                    if respFlagFix
                        for t = 1:size(metric.respFix,4)
                            respTmp = metric.respFix(:,:,:,t,cond);
                            respTmp(fitMask) = respFix(:,:,:,t,cond);
                            metric.respFix(:,:,:,t,cond) = respTmp;
                        end
                        clear respTmp
                    end
                    if respFlagMix
                        for t = 1:size(metric.respFix,4)
                            for run = 1:size(metric.respMix,6)
                                respTmp = metric.respMix(:,:,:,t,cond,run);
                                respTmp(fitMask) = respMix(:,:,:,t,run,cond);
                                metric.respMix(:,:,:,t,cond,run) = respTmp;
                            end
                        end
                    end
                    if respFlag || respFlag2
                        %temporal responses
                        for run = 1:nRuns
                            curRun = (cond-1)*nRuns+run;
                            if respFlag
                                tmpResp = mean(reshape(dataDetrend{curRun},[sz(1:3) 12 sz(4)/12]),5);
                                for t = 1:size(tmpResp,4)
                                    tmp = nan(size(metric.mask));
                                    tmp(ind) = tmpResp(:,:,:,t);
                                    metric.resp(:,:,:,t,cond,run) = tmp;
                                end
                            end
                            if respFlag2
                                tmpResp2 = dataDetrend{curRun};
                                for t = 1:size(tmpResp2,4)
                                    tmp = nan(size(metric.mask));
                                    tmp(ind) = tmpResp2(:,:,:,t);
                                    metric.resp2(:,:,:,t,cond,run) = tmp;
                                end
                            end
                        end
                    end
                end
                %F
                tmpF = fixedF;
                tmpFDR = fixedFDR;
                tmp = nan(size(metric.mask));
                tmp(ind) = tmpF;
                metric.allExp.F = tmp;
                tmp = nan(size(metric.mask));
                tmp(ind) = tmpFDR;
                metric.allExp.FDR = tmp;
                
                
                %Preselect session
                sessionLabelOrig = cell2mat(results.inputs.opt.sessionLabel);
                sessionLabel = sessionLabelOrig(1:end/3);
                switch p.preSelection
                    case {'sess1' 'sess2'}
                        metric.allRuns.amp = metric.allRuns.amp(:,:,:,sessionLabel==str2num(p.preSelection(5)),:);
                        metric.allRuns.delay = metric.allRuns.delay(:,:,:,sessionLabel==str2num(p.preSelection(5)),:);
                        if respFlag
                            metric.resp = metric.resp(:,:,:,:,:,sessionLabel==str2num(p.preSelection(5)));
                        end
                        if respFlag2
                            metric.resp2 = metric.resp2(:,:,:,:,:,sessionLabel==str2num(p.preSelection(5)));
                        end
                        nRuns = size(metric.allRuns.amp,4);
                        sessionLabel = ones(1,size(metric.allRuns.amp,4));
                    case 'bothSess'
                        %                         sessionLabel = ones(1,size(metric.allRuns.amp,4));
                        %                         sessionLabel(end/2+1:end) = 2;
                    otherwise
                        error('did not program this yet')
                        sessionLabel = ones(1,size(metric.allRuns.amp,4)); %%%%%%%%%%
                end % for metric.allExp, preSelection already done
                period = results.inputs.stimdur*2;
                clear results
                
                %Apply anatomical mask
                allFields = fields(metric); allFields(ismember(allFields, {'mask' 'brain','ecc','pol'})) = [];
                for i = 1:length(allFields)
                    if isstruct(metric.(allFields{i}))
                        allFields2 = fields(metric.(allFields{i}));
                        for ii = 1:length(allFields2)
                            %                         [allFields{i} '.' allFields2{ii}]
                            %                         size(metric.(allFields{i}).(allFields2{ii}))
                            %                         size(repmat(~metric.mask,[1 1 1 size(metric.(allFields{i}).(allFields2{ii}),4)]))
                            metric.(allFields{i}).(allFields2{ii})(repmat(~metric.mask,[1 1 1 size(metric.(allFields{i}).(allFields2{ii}),4)])) = nan;
                            %                         metric.(allFields{i}).(allFields2{ii})(repmat(~metric.mask,[1 1 1 size(metric.(allFields{i}).(allFields2{ii}),4),size(metric.(allFields{i}).(allFields2{ii}),5),size(metric.(allFields{i}).(allFields2{ii}),6)])) = nan;
                        end
                    else
                        switch allFields{i}
                            case {'resp' 'resp2','respFix','respMix'}
                                sz = size(metric.(allFields{i}));
                                metric.(allFields{i})(repmat(~metric.mask,[1 1 1 sz(4:end)])) = nan;
                            otherwise
                                error('XX')
                        end
                    end
                end
                
                %Define funcROI
                p.funcROI.stats = metric.allExp.F;
                
                
                
%                 %%%%%%%
%                 % create some image for publication
%                 close all
%                 f1 = figure('WindowStyle','docked');
%                 f2 = figure('WindowStyle','docked');
% %                 maxF = 20;
%                 maxF = 2*pi;
%                 for slice = 12%1:size(metric.brain,3)
%                     brain = metric.brain(:,:,slice);
% %                     overlay = metric.allExp.F(:,:,slice);
%                     overlay = metric.allRuns.delay(:,:,slice,1,1);
%                     overlay = wrapToPi(overlay-circ_mean(overlay(~isnan(overlay))))+pi;
%                     brain = imadjust(brain./max(max(brain)));
%                     figure(f1)
% %                     subplot(3,ceil(size(metric.brain,3)/3),slice)
%                     imagesc(overlay,[0 maxF]); colormap jet
%                     colorbar
%                     
%                     overlayMask = repmat(~isnan(overlay),[1 1 3]);
%                     
%                     brain = repmat(brain,[1 1 3]);
%                     brain = brain./max(max(max(brain)));
%                     brain = uint8(floor(brain * 255));
%                     overlay = overlay./maxF;
%                     overlay = uint8(floor(overlay* 255));
%                     overlay = ind2rgb(overlay, jet(255));
%                     overlay = uint8(floor(overlay* 255));
%                     brainOverlay = brain;
%                     brainOverlay(overlayMask) = overlay(overlayMask);
%                     figure(f2)
% %                     subplot(3,ceil(size(metric.brain,3)/3),slice)
%                     imshow(brainOverlay)
%                 end
%                 %%%%%%%
                
                
                
                %Make this fit the old pipeline
                metricSinCos = metric;
                metric = struct;
                
                metric.voxInd = [];
                metric.ampRect = nan(length(find(metricSinCos.mask)),nRuns,3);
                metric.delayRect = nan(length(find(metricSinCos.mask)),nRuns,3);
                metric.ecc = nan(length(find(metricSinCos.mask)),1);
                metric.pol = nan(length(find(metricSinCos.mask)),1);
                if respFlagFix
                    metric.respFix = nan(length(find(metricSinCos.mask)),12,3);
                end
                if respFlagMix
                    metric.respMix = nan(length(find(metricSinCos.mask)),12,3,nRuns);
                end
                if respFlag
                    metric.resp = nan(length(find(metricSinCos.mask)),12,nRuns,3);
                end
                if respFlag2
                    metric.resp2 = nan(length(find(metricSinCos.mask)),12*8,nRuns,3);
                end
                metric.funcROI.stats = nan(length(find(metricSinCos.mask)),1);
                for z = 1:size(metricSinCos.mask,3)
                    [x,y] = find(metricSinCos.mask(:,:,z));
                    metric.voxInd = [metric.voxInd; x y repmat(z,[length(x) 1])];
                end
                
                
                for ii = 1:size(metric.voxInd,1)
                    metric.ampRect(ii,:,:) = metricSinCos.allRuns.amp(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3),:,:);
                    metric.delayRect(ii,:,:) = metricSinCos.allRuns.delay(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3),:,:);
                    if respFlagFix
                        metric.respFix(ii,:,:) = permute(metricSinCos.respFix(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3),:,:),[1 4 5 2 3]);
                    end
                    if respFlagMix
                        metric.respMix(ii,:,:,:) = permute(metricSinCos.respMix(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3),:,:,:),[1 4 5 6 2 3]);
                    end
                    if respFlag
                        metric.resp(ii,:,:,:) = permute(metricSinCos.resp(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3),:,:,:),[1 4 6 5 2 3]);
                    end
                    if respFlag2
                        metric.resp2(ii,:,:,:) = permute(metricSinCos.resp2(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3),:,:,:),[1 4 6 5 2 3]);
                    end
                    metric.funcROI.stats(ii,1) = p.funcROI.stats(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3));
                    metric.ecc(ii,1) = metricSinCos.ecc(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3));
                    metric.pol(ii,1) = metricSinCos.pol(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3));
                    
                    metric.fixedAmpRect(ii,1,:) = metricSinCos.allExp.amp(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3),:,:);
                    metric.fixedDelayRect(ii,1,:) = metricSinCos.allExp.delay(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3),:,:);
                    metric.fixedF(ii,1,1) = metricSinCos.allExp.F(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3));
                    metric.fixedFDR(ii,1,1) = metricSinCos.allExp.FDR(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3));
                end
                
                %Sort functional ROI and put in p struct
                p.funcROI.vec.stats = metric.funcROI.stats;
                [~, p.funcROI.vec.sortInd] = sort(metric.funcROI.stats,'descend');
                
                %Define other stuff
                metric.sessionLabel(1,:,:) = repmat(sessionLabel',[1 3]);
                metric.fittingParam.period = period;
                metric.mask = metricSinCos.mask;
                metric.brain = metricSinCos.brain;
                
                
                clear metricSinCos
            otherwise
                error('did not implement for these specific params')
        end
        
        
        %% Define some params
        p.subj = subj;
        p.fitType = fitType;
        %         p.maskName = maskName;
        p.mask = metric.mask;
        p.fitMask = fitMask;
        p.brain = metric.brain;
        p.k = k;
        p.nFeatLev = nFeatLev;
        p.keepInfo = keepInfo;
        %         p.C = 1;
        p.voxInd = metric.voxInd;
        
        p.nObs = size(metric.ampRect,2)*2;
        p.nFeatures = size(metric.ampRect,1);
        
        p.dataFileOut = dataFileOut;
        p.dataDirOut = dataDirOut;
        
        switch p.featLevType
            case 'max'
                p.featLevList = p.nFeatures;
            case 'log'
                p.featLevList = unique([1 round(exp(linspace(log(1),log(p.nFeatures),p.nFeatLev))) p.nFeatures]);
            case 'linReg'
                p.featLevList = p.voxPerFeatLev:p.voxPerFeatLev:p.nFeatures;
                if p.nFeatures~=p.featLevList(end)
                    p.featLevList(end+1) = p.nFeatures;
                end
            case 'lin'
                error('doubleCheck')
                p.featLevList = linspace(0,p.nFeatures,p.nFeatLev+1);
                p.featLevList(1) = [];
                p.featLevList = round(p.featLevList);
            case 'F1'
                Fsorted = p.funcROI.vec.stats(p.funcROI.vec.sortInd);
                p.funLevList = 100;
                FatStep = Fsorted(p.funLevList);
                Fsorted = Fsorted(p.funLevList:end);
                Fstep = linspace(max(Fsorted),min(Fsorted),p.nFeatLev);
                for iii = 2:length(Fstep)
                    p.funLevList = [p.funLevList p.funLevList(1)+length(find(Fsorted>Fstep(iii)))];
                    FatStep = [FatStep Fstep(iii)];
                end
                p.featLevList = p.funLevList;
                
            case 'fix'
                p.featLevList = [1 2 3 4 5 6 7 8 9 10 12 14 16 19 22 25 29 33 39 45 52 60 69 80 93 108 125 144 167 193 224 259 300 347 401 465 538 622 720 834 965 1117 1293];
                p.featLevList = p.featLevList(p.featLevList<=p.nFeatures);
                if p.featLevList(end)~=p.nFeatures
                    p.featLevList = [p.featLevList p.nFeatures];
                end
            case 'fix2'
                p.featLevList = [10 100 500 1000];
                p.featLevList = p.featLevList(p.featLevList<=p.nFeatures);
                if p.featLevList(end)~=p.nFeatures
                    p.featLevList = [p.featLevList p.nFeatures];
                end
            otherwise
                error('what do you wanna do?')
        end
        p.nFeatLev = length(p.featLevList);
        
        
        
        %         if strfind(char(curKernel),'polynomial')
        %             tmpKernel = strsplit(char(curKernel),'_');
        %             p.kernel = tmpKernel{1}; % 'linear' 'quadratic' 'polynomial' 'rbf' (Gaussian Radial) 'mlp' (Multilayer Perceptron)
        %             p.polyorder = str2num(tmpKernel{2});
        %         else
        %             p.kernel = char(curKernel); % 'linear' 'quadratic' 'polynomial' 'rbf' (Gaussian Radial) 'mlp' (Multilayer Perceptron)
        %         end
        
        %% Define data and labels
        %Data
        switch p.conditions
            case 'grat1VSgrat2'
                ind1 = 1; ind2 = 2; indNorm = 3;
            case 'gratVSplaid'
                ind1 = [1 2]; ind2 = 3; indNorm = 3;
            case 'grat1VSplaid'
                ind1 = 1; ind2 = 3; indNorm = 2;
            case 'grat2VSplaid'
                ind1 = 2; ind2 = 3; indNorm = 1;
            otherwise
                error('something wrong')
        end
        if strcmp(fitType,'SinCos')
            if ~strcmp(p.conditions,'gratVSplaid')
                %mixed-effect
                tmp1 = metric.ampRect(:,:,ind1);
                tmp2 = metric.ampRect(:,:,ind2);
                tmpNorm = metric.ampRect(:,:,indNorm);
                mag = [tmp1 tmp2]';
                magNorm = tmpNorm';
                
                tmp1 = metric.delayRect(:,:,ind1);
                tmp2 = metric.delayRect(:,:,ind2);
                tmpNorm = metric.delayRect(:,:,indNorm);
                phase = [tmp1 tmp2]';
                phaseNorm = tmpNorm';
                
                %fixed-effect
                tmp1 = metric.fixedAmpRect(:,:,ind1);
                tmp2 = metric.fixedAmpRect(:,:,ind2);
                tmpNorm = metric.fixedAmpRect(:,:,indNorm);
                magFixed = [tmp1 tmp2]';
                magNormFixed = tmpNorm';
                
                tmp1 = metric.fixedDelayRect(:,:,ind1);
                tmp2 = metric.fixedDelayRect(:,:,ind2);
                tmpNorm = metric.fixedDelayRect(:,:,indNorm);
                phaseFixed = [tmp1 tmp2]';
                phaseNormFixed = tmpNorm';
            else
                %Average on polar data
                tmp1 = metric.ampRect(:,:,ind1);
                tmp2 = metric.ampRect(:,:,ind2);
                tmpNorm = metric.ampRect(:,:,indNorm);
                mag = [mean(tmp1,3) tmp2]';
                magNorm = tmpNorm';
                
                tmp1 = metric.delayRect(:,:,ind1);
                tmp2 = metric.delayRect(:,:,ind2);
                tmpNorm = metric.delayRect(:,:,indNorm);
                phase = [circ_mean(tmp1,[],3) tmp2]';
                phaseNorm = tmpNorm';
                %                 %Averge on cartesian data
                %                 [X,Y] = pol2cart(metric.delayRect,metric.ampRect);
                %                 tmp1 = mean(X(:,:,ind1),3);
                %                 tmp2 = X(:,:,ind2);
                %                 tmpNorm = X(:,:,indNorm);
                %                 X = [tmp1 tmp2]';
                %                 XNorm = tmpNorm';
                %
                %                 tmp1 = mean(Y(:,:,ind1),3);
                %                 tmp2 = Y(:,:,ind2);
                %                 tmpNorm = metric.delayRect(:,:,indNorm);
                %                 Y = [tmp1 tmp2]';
                %                 YNorm = tmpNorm';
                %
                %                 [phase,mag] = cart2pol(X,Y);
                %                 [phaseNorm,magNorm] = cart2pol(XNorm,YNorm);
            end
        else
            error('have a look at that (norm)')
            tmp1 = mean(metric.ampRect(:,:,ind1),3);
            tmp2 = mean(metric.ampRect(:,:,ind2),3);
            mag = [tmp1 tmp2]';
            tmp1 = circ_mean(sec2rad(metric.delayRect(:,:,ind1),metric.fittingParam.period),[],3);
            tmp2 = circ_mean(sec2rad(metric.delayRect(:,:,ind2),metric.fittingParam.period),[],3);
            phase = [tmp1 tmp2]';
        end
        [X,Y] = pol2cart(phase,mag);
        d.xData = complex(X,Y);
        [X,Y] = pol2cart(phaseNorm,magNorm);
        d.normData = complex(X,Y); clear tmp1 tmp2 tmpNorm phase phaseNorm mag magNorm
        [X,Y] = pol2cart(phaseFixed,magFixed);
        d.xDataFixed = complex(X,Y);
        [X,Y] = pol2cart(phaseNormFixed,magNormFixed);
        d.normDataFixed = complex(X,Y); clear phaseFixed phaseNormFixed magFixed magNormFixed
        md.voxInd = metric.voxInd';
        md.mask = metric.mask;
        md.brain = metric.brain;
        
        %Category labels
        d.label = reshape(repmat([1 2],p.nObs/2,1),p.nObs,1);
        d.labelFixed = reshape(repmat([1 2],1,1),2,1);
        %Session labels
        if length(ind1)>1 && ~all(metric.sessionLabel(1,:,ind1(1))==metric.sessionLabel(1,:,ind1(2)))
            error('something wrong here')
        end
        if length(ind2)>1
            error('something wrong here')
        end
        d.sessionLabel = [metric.sessionLabel(1,:,ind1(1)) metric.sessionLabel(1,:,ind2)]';
        d.sessionLabelFixed = [metric.sessionLabel(1,1,ind1(1)) metric.sessionLabel(1,1,ind2)]';
        
        %F
        d.fixedF = metric.fixedF;
        d.fixedFDR = metric.fixedFDR;
        
        %ecc and pol
        d.ecc = metric.ecc;
        d.pol = metric.pol;
        
        %% Compute accuracy threshold on binomial distribution !!!!Might be wrong!!!!
        [p.binom.hitRate, p.binom.p] = binomialSigThresh(p.nObs,0.5);
        
        
        
        
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %% Test of independence (do the thing on half the data)
        %         startAt = 2;
        %         d.xData = d.xData(startAt:2:end,:);
        %         d.normData = d.normData(startAt:2:end,:);
        %         d.label = d.label(startAt:2:end,:);
        %         d.sessionLabel = d.sessionLabel(startAt:2:end,:);
        %         p.nObs = p.nObs/2;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        p.nTrialPerRun = str2double(p.split(1));
        if p.splitHalf
            error('double-check with respect to fixed')
            d.xData = d.xData(p.splitHalf:2:end,:);
            d.normData = d.normData(p.splitHalf:2:end,:);
            d.label = d.label(p.splitHalf:2:end,:);
            d.sessionLabel = d.sessionLabel(p.splitHalf:2:end,:);
            p.nObs = p.nObs/2;
            p.nTrialPerRun = p.nTrialPerRun/2;
        end
        d.runLabel = reshape(repmat(1:p.nObs /p.nTrialPerRun,p.nTrialPerRun,1),p.nObs,1);
        
        %         if length(p.infoComb)==1
        %             error('XX')
        %             switch p.infoComb{1}
        %                 case {'catComb','infoComb_Copt','voxComb_Copt'}
        %                 case 'pcaOnTS_lowLevel'
        %                     sz = size(metric.resp2(:,:,:,1:2));
        %                     d.xData = permute(reshape(metric.resp2(:,:,:,1:2),[sz(1:2) prod(sz(3:4))]),[3 1 2]);
        %                     d.normData = [];
        %                 case {'pcaOnResp_lowLevel','pcaOnResp_pcaOnSpace'}
        %                     sz = size(metric.resp(:,:,:,1:2));
        %                     d.xData = permute(reshape(metric.resp(:,:,:,1:2),[sz(1:2) prod(sz(3:4))]),[3 1 2]);
        %                     d.normData = [];
        %
        %                 case {'respSVM__pcaNone_pcaNone','respSVM__pcaRun_pcaNone','respSVM__pcaTime_pcaNone','respSVM__pcaNone_pcaSpace','respSVM__pcaRun_pcaSpace','respSVM__pcaTime_pcaSpace','respSVM__pcaNone_pcaRun','respSVM__pcaRun_pcaRun','respSVM__pcaTime_pcaRun'}
        %                     sz = size(metric.resp(:,:,:,1:2)); % vox x time x run x cond
        %                     d.xData = permute(reshape(metric.resp(:,:,:,1:2),[sz(1:2) prod(sz(3:4))]),[3 1 2]);  % run x vox x time
        %                     d.normData = [];
        %
        %                 otherwise
        %                     error('XX')
        %             end
        %         else
        if any(ismember(p.infoComb,tmpStrResp))
            sz = size(metric.resp(:,:,:,1:2)); % vox x time x run x cond
            d.respData = permute(reshape(metric.resp(:,:,:,1:2),[sz(1:2) prod(sz(3:4))]),[3 1 2]);  % run x vox x time
            [sz(1), sz(2), sz(3), sz(4)] = size(metric.resp(:,:,:,3)); % vox x time x run x cond
            d.respDataNorm = permute(reshape(metric.resp(:,:,:,3),[sz(1:2) prod(sz(3:4))]),[3 1 2]);  % run x vox x time
        else
            warning('resp not loaded')
        end
        %         end
        
        if isfield(p,'saveDataOnly') && p.saveDataOnly
%             save(fullfile(dataDirOut,[dataFileOut '__data.mat']),'md','d','p','Fmask')
            if respFlagFix
                resp.data.fix = metric.respFix;
            end
            if respFlagMix
                resp.data.mix = metric.respMix;
            end
            resp.info = 'vox(as in d struct) x t x cond x run';
            save(fullfile(dataDirOut,[p.subj '_' p.preSelection '.mat']),'d','resp','p')
            save(fullfile(dataDirOut,[p.subj '_' p.preSelection '_metric.mat']),'metric')
            continue
        end
        
        
        
        
        %% PCA denoise
        if strcmp(p.infoComb,'respPCA')
%             %Run spatial PCA on timeseries
%             sz2= size(metric.resp2);
%             ts = reshape(metric.resp2,[sz2(1:2) prod(sz2(3:4))]);
%             sz2= size(ts);
%             ts2 = reshape(ts,[sz2(1) prod(sz2(2:3))])';
%             [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(ts2);
%             %Fit PCs with sine and cosine
%             %will be more easily done earlier...
        end

        
        
        %% Run SVMs
        p.doPerm = 0;
        stime = tic;
        if strcmp(p.k,'none')
            error('double-check the inclusion of ampAtDelay and ampAtVoxDelay')
            if perm.doPerm
                svm = runSVMnoK(d,p,1,ppn);
            else
                svm = runSVMnoK(d,p,1,0);
            end
        else
%                         svm = runSVMRepetitions4(d,p,1,ppn);
            svm = runSVMRepetitions5(d,p,1,ppn);
            if isempty(svm)
                close all
                continue
            end
        end
        
        
        
        %% Output and plot
        if stopThis
            break
        else
            %Save
            try
                save(fullfile(dataDirOut,[dataFileOut '.mat']),'svm','md','Fmask')
            catch
                save(fullfile(dataDirOut,[dataFileOut '.mat']),'svm','md','Fmask','-v7.3')
            end
            display(['Data saved as ' dataFileOut])
            
            if p.plotIt
                %Plot
                %             h = plotSVM_forSVMOneRepetition2(svm);
                %             h = plotSVM_forSVMOneRepetition3(svm);
                if length(p.infoComb)==1 && strcmp(p.infoComb,'crossValDelay')
                    h = plotHistAndDelay(svm.rTr.diff,svm.rTe.diff);
                    xlabel('delay difference (rad)')
                    title([svm.p.subj '; ' svm.p.preSelection])
                    saveas(h,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{1} '_diff']),'jpg')
                    saveas(h,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{1} '_diff']))
                    %                     h = plotHistAndDelay(svm.rTr.t,svm.rTe.t);
                    %                     xlabel('t-value for delay difference')
                    %                     saveas(h,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{1} '_tDiff']),'jpg')
                    %                     saveas(h,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{1} '_tDiff']))
                elseif length(p.infoComb)==1 && strcmp(p.infoComb,'crossValDelayQ')
                    edges = circ_mean(svm.rTr.diff.edges,[],1);
                    binWidth = diff(edges);
                    binCent = edges(1:end-1)+binWidth/2;
                    N = mean(svm.rTr.diff.N,1);
                    delay = circ_mean(svm.rTr.diff.delay,[],1);
                    crossValDelay = circ_mean(svm.rTe.diff.delay,[],1);
                    quant = linspace(0+100/length(binCent)/2,100-100/length(binCent)/2,length(binCent));

%                     figure('WindowStyle','docked');
%                     yyaxis right
%                     plot(binCent./pi*6,delay./pi*6,'-o'); hold on
%                     plot(binCent./pi*6,crossValDelay./pi*6,'-*');
%                     plot(get(gca,'xlim'),[0 0],':')
%                     yyaxis left
%                     bar(binCent./pi*6,N);
                    
                    figure('WindowStyle','docked');
                    yyaxis right
                    h(1) = plot(quant,delay./pi*6,'-o'); hold on
                    h(2) = plot(quant,crossValDelay./pi*6,'-*');
                    ylabel('delay (sec)')
                    plot(get(gca,'xlim'),[0 0],':')
                    yyaxis left
                    bar(quant,N);
                    ylabel('vox count')
                    legend(h,{'delay' 'cross-calidated delay'})
                    xlabel('quantile (%)')

                else
                    
                    if isfield(svm,'rTr') && isfield(svm.rTr{1},'cOpt')
                        hCC = figure('WindowStyle','docked');
                        for anaIndX = 1:length(svm.rTr{1}.cOpt)
                            subplot(1,length(svm.rTr{1}.cOpt),anaIndX)
                            [C,hC] = contour(log(svm.p.cList2),log(svm.p.cList),svm.rTr{1}.cOpt{anaIndX}.distT(:,:,1));
                            clabel(C,hC)
                            if anaIndX==1
                                ylabel('high-level C')
                                xlabel('low-level C')
                            end
                            title(svm.rTr{1}.cOpt{anaIndX}.info)
                        end
                        saveas(hCC,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__contour.jpg']),'jpg')
                        saveas(hCC,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__contour.fig']),'fig')
                        close(hCC)
                    end
                    
                    
                    if length(svm.rTr)==1 && isfield(svm.rTr{1},'latent')
                        plotPCA(svm.rTr{1},svm.p)
                        figure('WindowStyle','docked');
                        distT = mean(svm.rTe{1}.distT,4);
                        distT_std = std(svm.rTe{1}.distT,[],4);
                        bar(distT); hold on
%                         errorbar(distT,distT_std,'.b');
                        set(gca,'Xtick',1:length(distT))
                        xlabel('spatial PCs')
                        ylabel('distT')
                        title([svm.p.subj '; ' svm.p.preSelection])
                        saveas(gca,fullfile(p.dataDirOut,[p.dataFileOut '__spatialPCA_pred.jpg']))
                    else
                        h = plotSVM_forParaAcc(svm,[],[],0);
                        if length(h)==1
                            saveas(h,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '.jpg']),'jpg')
                            close(h)
                        else
                            for anaInd = 1:length(h)
                                saveas(h(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '.jpg']),'jpg')
                                close(h(anaInd))
                            end
                        end
                        
                        h2 = plotSVM_forParaAcc(svm,[],[],1);
                        if length(h2)==1
                            saveas(h2,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '_abs.jpg']),'jpg')
                            close(h2)
                        else
                            for anaInd = 1:length(h2)
                                saveas(h2(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_abs.jpg']),'jpg')
                                close(h2(anaInd))
                            end
                        end
                        
                        h3 = plotSVM_forParaDist(svm,[],[],1);
                        if length(h3)==1
                            saveas(h3,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '_distAbs.jpg']),'jpg')
                            close(h3)
                        else
                            for anaInd = 1:length(h3)
                                saveas(h3(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_distAbs.jpg']),'jpg')
                                close(h3(anaInd))
                            end
                        end
                        
                        h4 = plotSVM_forParaDist(svm,[],[],1,1);
                        if length(h4)==1
                            saveas(h4,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '_distTabs.jpg']),'jpg')
                            close(h4)
                        else
                            for anaInd = 1:length(h4)
                                saveas(h4(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_distTabs.jpg']),'jpg')
                                close(h4(anaInd))
                            end
                        end
                        
                        h5 = plotSVM_forParaDist(svm,[],[],0);
                        if length(h5)==1
                            saveas(h5,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '_dist.jpg']),'jpg')
                            close(h5)
                        else
                            for anaInd = 1:length(h5)
                                saveas(h5(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_dist.jpg']),'jpg')
                                close(h5(anaInd))
                            end
                        end
                        
                        h6 = plotSVM_forParaDist(svm,[],[],0,1);
                        if length(h6)==1
                            saveas(h6,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '_distT.jpg']),'jpg')
                            close(h6)
                        else
                            for anaInd = 1:length(h6)
                                saveas(h6(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_distT.jpg']),'jpg')
                                close(h6(anaInd))
                            end
                        end
                    end
                end
            end
        end
        
        
        svm.dur = toc(stime);
        display(['This took ' datestr(1/24/60/60*svm.dur, 'HH:MM:SS') 'sec to compute, compile rep and save'])
        
    else
        %% Load already performed SVM analysis
        if length(svm.rTr)==1 && isfield(svm.rTr{1},'latent')
%             plotPCA(svm.rTr{1},svm.p)
            h = figure('WindowStyle','docked');
            distT = mean(svm.rTe{1}.distT,4);
            distT_std = std(svm.rTe{1}.distT,[],4);
            bar(distT); hold on
            set(gca,'Xtick',1:length(distT))
            xlabel('spatial PCs')
            ylabel('distT')
            if ~isfield(svm,'perm')
                title([svm.p.subj '; ' svm.p.preSelection])
            else
                title([svm.p.subj '; ' svm.p.preSelection '; ' num2str(size(svm.perm.rTe{1}.distT,4)) 'perms'])
                distT2 = mean(svm.perm.rTe{1}.distT,4);
                distT2_std = squeeze(prctile(svm.perm.rTe{1}.distT,[5 95],4));
                errorbar(1:length(distT2),distT2,distT2_std(:,1)'-distT2,distT2_std(:,2)'-distT2,'or');
            end
        else
            load(allFile{taskInd},'svm','md','Fmask')
            h = plotSVM_forParaAcc(svm,[],[],0);
%             h2 = plotSVM_forParaAcc(svm,[],[],1);
%             h3 = plotSVM_forParaDist(svm,[],[],1);
%             h4 = plotSVM_forParaDist(svm,[],[],1,1);
%             h5 = plotSVM_forParaDist(svm,[],[],0);
%             h6 = plotSVM_forParaDist(svm,[],[],0,1);
        end

%         h = plotSVM_forParaDis(svm);
%         h = plotSVM_forSVMOneRepetition2(svm);
%         h = plotSVM_forSVMOneRepetition3(svm); % for all in different subplots
        drawnow;
    end
    
    %% Permutation test
    if ~isempty(perm) && perm.doPerm
        d = svm.d;
        p = svm.p;
        p.perm = perm;
        p.doPerm = 1;
        p.repeat = perm.compileEachNperm;
        svm.p = p;
        if ~isfield(svm,'perm')
            svm.perm = svm;
            svm.perm.rTr = clearAllDataField2(svm.perm.rTr,[],1);
            svm.perm.rVal = clearAllDataField2(svm.perm.rVal,[],1);
            svm.perm.rTe = clearAllDataField2(svm.perm.rTe,[],1);
        end
%         svmCompiled = svm;
%         svm = svmCompiled.perm;
        
        
        % Loop over permutations in the bloc
        if ~isfield(svm.perm.rTe{1},'acc')
            if size(svm.perm.rTe{1}.distT,4)<p.perm.numPerm
                permDone = false;
            else
                permDone = true;
            end
        else
            if length(svm.perm.rTe{1}.acc)<p.perm.numPerm
                permDone = false;
            else
                permDone = true;
                %             %Plot
                %             h = plotSVM_forSVMOneRepetition3(svm,h); % for all in different plots
            end
        end
        
        i=0;
        while i<svm.p.perm.numPermPerBloc && ~permDone
            svm_cur = runSVMRepetitions5(d,p,1,ppn);

            tic
            display(['Compiling and saving as ' svm.p.dataFileOut])
            % Compile permutations in bloc
            svm_cur.rTr = clearAllDataField2(svm_cur.rTr,[],0);
            svm_cur.rVal = clearAllDataField2(svm_cur.rVal,[],0);
            svm_cur.rTe = clearAllDataField2(svm_cur.rTe,[],0);
            
            svm.perm.rTr = compileRep2(svm.perm.rTr,svm_cur.rTr);
            svm.perm.rVal = compileRep2(svm.perm.rVal,svm_cur.rVal);
            svm.perm.rTe = compileRep2(svm.perm.rTe,svm_cur.rTe);
            
            
            if ~isfield(svm.perm.rTe{1},'acc')
                sz = size(svm_cur.rTe{1}.distT);
                i = i+size(svm_cur.rTe{1}.distT,length(sz));
                if size(svm_cur.rTe{1}.distT,length(sz))>=p.perm.numPerm
                    permDone = true;
                end
            else
                sz = size(svm_cur.rTe{1}.acc);
                i = i+size(svm_cur.rTe{1}.acc,length(sz));
                if size(svm.perm.rTe{1}.acc,length(sz))>=p.perm.numPerm
                    permDone = true;
                end
            end
            
            
            if length(svm.rTr)==1 && isfield(svm.rTr{1},'latent')
                hold off
                bar(distT); hold on
                distT2 = mean(svm.perm.rTe{1}.distT,4);
                distT2_std = squeeze(prctile(svm.perm.rTe{1}.distT,[5 95],4));
                errorbar(1:length(distT2),distT2,distT2_std(:,1)'-distT2,distT2_std(:,2)'-distT2,'or');
                set(gca,'Xtick',1:length(distT))
                xlabel('spatial PCs')
                ylabel('distT')
                title([svm.p.subj '; ' svm.p.preSelection '; ' num2str(size(svm.perm.rTe{1}.distT,4)) 'perms'])
                saveas(gca,fullfile(p.dataDirOut,[p.dataFileOut '__spatialPCA_pred.jpg']))
            else
                %Plot
                %             h = plotSVM_forSVMOneRepetition2(svm,h);
                %             h = plotSVM_forSVMOneRepetition3(svm,h); % for all in different plots
                h = plotSVM_forParaAcc(svm,h,[],0);
                for anaInd = 1:length(h)
                    saveas(h(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '.jpg']),'jpg')
                end
%                 h2 = plotSVM_forParaAcc(svm,h2,[],1);
%                 for anaInd = 1:length(h2)
%                     saveas(h2(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_abs.jpg']),'jpg')
%                 end
%                 h3 = plotSVM_forParaDist(svm,h3,[],1);
%                 for anaInd = 1:length(h3)
%                     saveas(h3(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_distAbs.jpg']),'jpg')
%                 end
%                 h4 = plotSVM_forParaDist(svm,h4,[],1,1);
%                 for anaInd = 1:length(h4)
%                     saveas(h4(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_distTabs.jpg']),'jpg')
%                 end
%                 h5 = plotSVM_forParaDist(svm,h5,[],0);
%                 for anaInd = 1:length(h5)
%                     saveas(h5(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_dist.jpg']),'jpg')
%                 end
%                 h6 = plotSVM_forParaDist(svm,h6,[],0,1);
%                 for anaInd = 1:length(h6)
%                     saveas(h6(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_distT.jpg']),'jpg')
%                 end
            end

            %Save
            try
                save(fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '.mat']),'svm','-append')
            catch
                save(fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '.mat']),'svm','-append','-v7.3')
            end
            display(['This took ' datestr(1/24/60/60*toc, 'HH:MM:SS') 'sec to compute, compile rep and save'])
            
        end
        for anaInd = 1:length(h)
            if exist('h1','var'); close(h(anaInd)); end
            if exist('h2','var'); close(h2(anaInd)); end
            if exist('h3','var'); close(h3(anaInd)); end
            if exist('h4','var'); close(h4(anaInd)); end
            if exist('h5','var'); close(h5(anaInd)); end
            if exist('h6','var'); close(h6(anaInd)); end
        end
        
        taskDone(taskInd) = permDone;
    end
end





function [FDR,pGood] = getPfromF(F,df)
%Convert F to p
%Convert to p and Z, with trick to avoid p=0 and z=+inf at very large F
p = fcdf(F,df.pFull-df.pReduced,df.n-df.pFull);
q = fcdf(1./F,df.n-df.pFull,df.pFull-df.pReduced);
pGood = nan(size(p));
pGood(p<=0.5) = 1-p(p<=0.5);
pGood(p>0.5) = q(p>0.5);
%         [FDR,Q,PIO] = mafdr(pGood);
FDR = mafdr(pGood(:),'BHFDR',true);
tmp = nan(size(pGood));
tmp(1:numel(pGood)) = FDR;
FDR = tmp; clear tmp



function h = plotHistAndDelay(rTr,rTe)
binCent = rTr.binCent;
binCent2 = rTe.binCent2;
N = mean(rTr.N,1);
for i = 1:size(rTe.crossValDelay,2)
    if isempty(find(~isnan(rTe.crossValDelay(:,i)),1))
        crossValDelay(1,i) = nan;
        crossValDelay_std(1,i) = nan;
    else
        crossValDelay(1,i) = circ_mean(rTe.crossValDelay(~isnan(rTe.crossValDelay(:,i)),i),[],1);
        crossValDelay_std(1,i) = circ_std(rTe.crossValDelay(~isnan(rTe.crossValDelay(:,i)),i),[],[],1);
    end
end
for i = 1:size(rTe.crossValDelay2,2)
    if isempty(find(~isnan(rTe.crossValDelay2(:,i)),1))
        crossValDelay2(1,i) = nan;
        crossValDelay2_std(1,i) = nan;
    else
        crossValDelay2(1,i) = circ_mean(rTe.crossValDelay2(~isnan(rTe.crossValDelay2(:,i)),i),[],1);
        crossValDelay2_std(1,i) = circ_std(rTe.crossValDelay2(~isnan(rTe.crossValDelay2(:,i)),i),[],[],1);
    end
end


xLim = [rTr.binCent(1)-mean(diff(rTr.binCent))/2 rTr.binCent(end)+mean(diff(rTr.binCent))/2];
h = figure('WindowStyle','docked');
yyaxis left
bar(binCent,N,1,'faceColor',get(gca,'YColor'),'edgeColor','w'); hold on
ylabel('voxel count')
xlim(xLim)
yyaxis right
errorbar(binCent2(2:end-1),crossValDelay2(2:end-1)/pi*6,crossValDelay2_std(2:end-1)/pi*6,'-o'); hold on
errorbar(binCent2(1),crossValDelay2(1)/pi*6,crossValDelay2_std(1)/pi*6,'-o');
errorbar(binCent2(end),crossValDelay2(end)/pi*6,crossValDelay2_std(end)/pi*6,'-o');
ylabel('cross-validated delay difference (sec)')
plot(get(gca,'xlim'),[0 0],':')
