function [ output_args ] = runTask7(allFile,altPerm,altMainDir,guillimin,ppn)
rng('shuffle')
if ~exist('guillimin','var')
    guillimin=0;
end
if ~exist('ppn','var')
    ppn=0;
end

for taskInd = 1:length(allFile)
    clearvars -except taskInd allFile plotIt doPerm dataDirOut altPerm altMainDir guillimin ppn
    display(['Performing ' char(allFile{taskInd})])
    
    if length(find(size(allFile)>1))>1
        load(allFile{taskInd,1})
        load(allFile{taskInd,2})
    else
        load(allFile{taskInd})
    end
    
    %Load the good svm toolbox
    switch p.algorithm
        case {'runSVM_patAct2','runSVM_patAct4','runSVM_patAct5','runSVM_patAct6','runSVM_patAct7','runSVM_patAct8','runSVMOneRepetition','runSVMOneRepetition2'}
            addpath(genpath('C:\Users\Sebastien\OneDrive - McGill University\MATLAB\libsvm-3.21'))
        otherwise
            error('wanna use libsvm for that?')
            rmpath(genpath('C:\Users\Sebastien\OneDrive - McGill University\MATLAB\libsvm-3.21'))
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
                metric.mask = flipdim(permute(a.img,[3 1 2 4]),1); clear a
                metric.mask(:,:,1) = zeros(size(metric.mask,1),size(metric.mask,2));% Remove corrupted slices
                metric.mask(:,:,end) = zeros(size(metric.mask,1),size(metric.mask,2));% Remove corrupted slices
                a = load_nii(fullfile(dataDirIn,'run1a_preprocessing',ls(fullfile(dataDirIn,'run1a_preprocessing','trun101_preprocessed.nii*'))));
                metric.brain = mean(flipdim(permute(a.img,[3 1 2 4]),1),4); clear a
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Should review that if want to include ecc
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Apply F thresh to mask
                load(fullfile(dataDirIn,'sepSess170108_v1SinCos_1perRun_move12.mat'),'results')
                fitMask = results.mask;
                F = results.OLS.fixed.(p.preSelection).F.val.F;
                anatMask = logical(metric.mask);
                anatMask = anatMask(:,:,squeeze(any(any(fitMask,1),2)),:);
                anatMask = anatMask(:,squeeze(any(any(fitMask,1),3)),:,:);
                anatMask = anatMask(squeeze(any(any(fitMask,2),3)),:,:,:);
                
                [FDR,~] = getPfromF(F(logical(anatMask)),results.OLS.fixed.(p.preSelection).F.df);
                tmpFDR = nan(size(F));
                tmpFDR(logical(anatMask)) = FDR; FDR = tmpFDR;
                
                tmpFDR = nan(size(fitMask));
                tmpFDR(fitMask) = FDR; FDR = tmpFDR;
                
                Fmask = FDR>0.05;
                
                metric.mask(Fmask) = 0;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                if p.useEccMask
                    %Load ecc mask
                    a = load_nii(fullfile(dataDirIn,'run1c_masking',['lh.ecc.nii.gz']));
                    eccL = flipdim(permute(a.img,[3 1 2 4]),1); clear a
                    eccL(~metric.mask)=0;
                    a = load_nii(fullfile(dataDirIn,'run1c_masking',['rh.ecc.nii.gz']));
                    eccR = flipdim(permute(a.img,[3 1 2 4]),1); clear a
                    eccR(~metric.mask)=0;
                    eccMask = zeros(size(metric.mask));
                    eccMask(eccL>=p.eccLowLim&eccL<p.eccHighLim) = 1;
                    eccMask(eccR>=p.eccLowLim&eccR<p.eccHighLim) = 1;
                    metric.mask(~eccMask) = 0;
                end
                
                %Load data of interest
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 load(fullfile(dataDirIn,['161212_' upper(p.maskName) fitType '_' p.split '_' p.motionParam p.sm '.mat']),'results')
                %                 load(fullfile(dataDirIn,['161212_v1v2v3' fitType '_' p.split '_' p.motionParam p.sm '.mat']),'results')
                load(fullfile(dataDirIn,['X161212_v1' fitType '_' p.split '_' p.motionParam p.sm '.mat']),'results')
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
                    case 'sess2'
                        Ffield2 = 'Fcond3sess2';
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
                for cond = 1:3
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
                end
                
                %Preselect session
                sessionLabel = cell2mat(results.inputs.opt.sessionLabel);
                sessionLabel = sessionLabel(1:end/3);
                switch p.preSelection
                    case 'sess1'
                        metric.allRuns.amp = metric.allRuns.amp(:,:,:,sessionLabel==1,:);
                        metric.allRuns.delay = metric.allRuns.delay(:,:,:,sessionLabel==1,:);
                        nRuns = size(metric.allRuns.amp,4);
                        sessionLabel = ones(1,size(metric.allRuns.amp,4));
                    case 'sess2'
                        metric.allRuns.amp = metric.allRuns.amp(:,:,:,sessionLabel==2,:);
                        metric.allRuns.delay = metric.allRuns.delay(:,:,:,sessionLabel==2,:);
                        nRuns = size(metric.allRuns.amp,4);
                        sessionLabel = ones(1,size(metric.allRuns.amp,4));
                    case 'bothSess'
                        %                         sessionLabel = ones(1,size(metric.allRuns.amp,4));
                        %                         sessionLabel(end/2+1:end) = 2;
                    otherwise
                        error('did not program this yet')
                        sessionLabel = ones(1,size(metric.allRuns.amp,4)); %%%%%%%%%%
                end
                period = results.inputs.stimdur*2;
                clear results
                
                %Apply anatomical mask
                allFields = fields(metric); allFields(ismember(allFields, {'mask' 'brain'})) = [];
                for i = 1:length(allFields)
                    allFields2 = fields(metric.(allFields{i}));
                    for ii = 1:length(allFields2)
                        metric.(allFields{i}).(allFields2{ii})(repmat(~metric.mask,[1 1 1 size(metric.(allFields{i}).(allFields2{ii}),4)])) = nan;
                    end
                end
                
                %Define funcROI
                p.funcROI.stats = metric.allExp.F;
                
                %Make this fit the old pipeline
                metricSinCos = metric;
                metric = struct;
                
                metric.voxInd = [];
                metric.ampRect = nan(length(find(metricSinCos.mask)),nRuns,3);
                metric.delayRect = nan(length(find(metricSinCos.mask)),nRuns,3);
                metric.funcROI.stats = nan(length(find(metricSinCos.mask)),1);
                for z = 1:size(metricSinCos.mask,3)
                    [x,y] = find(metricSinCos.mask(:,:,z));
                    metric.voxInd = [metric.voxInd; x y repmat(z,[length(x) 1])];
                end
                
                
                for ii = 1:size(metric.voxInd,1)
                    metric.ampRect(ii,:,:) = metricSinCos.allRuns.amp(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3),:,:);
                    metric.delayRect(ii,:,:) = metricSinCos.allRuns.delay(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3),:,:);
                    metric.funcROI.stats(ii,1) = p.funcROI.stats(metric.voxInd(ii,1),metric.voxInd(ii,2),metric.voxInd(ii,3));
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
        md.voxInd = metric.voxInd';
        md.mask = metric.mask;
        md.brain = metric.brain;
        
        %Category labels
        d.label = reshape(repmat([1 2],p.nObs/2,1),p.nObs,1);
        %Session labels
        if length(ind1)>1 && ~all(metric.sessionLabel(1,:,ind1(1))==metric.sessionLabel(1,:,ind1(2)))
            error('something wrong here')
        end
        if length(ind2)>1
            error('something wrong here')
        end
        d.sessionLabel = [metric.sessionLabel(1,:,ind1(1)) metric.sessionLabel(1,:,ind2)]';
        
        %% Compute accuracy threshold on binomial distribution !!!!Might be wrong!!!!
        [p.binom.hitRate, p.binom.p] = binomialSigThresh(p.nObs,0.5);
        
        
        if isfield(p,'saveDataOnly') && p.saveDataOnly
            save(fullfile(dataDirOut,[dataFileOut '__data.mat']),'md','d','p','Fmask')
            continue
        end
        
        
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %% Test of independence (do the thing on half the data)
        %         startAt = 2;
        %         d.xData = d.xData(startAt:2:end,:);
        %         d.normData = d.normData(startAt:2:end,:);
        %         d.label = d.label(startAt:2:end,:);
        %         d.sessionLabel = d.sessionLabel(startAt:2:end,:);
        %         p.nObs = p.nObs/2;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if p.splitHalf
            d.xData = d.xData(p.splitHalf:2:end,:);
            d.normData = d.normData(p.splitHalf:2:end,:);
            d.label = d.label(p.splitHalf:2:end,:);
            d.sessionLabel = d.sessionLabel(p.splitHalf:2:end,:);
            p.nObs = p.nObs/2;
        end
        %% Run SVMs
        stime = tic;
        if strcmp(p.k,'none')
            error('double-check the inclusion of ampAtDelay and ampAtVoxDelay')
            if perm.doPerm
                svm = runSVMnoK(d,p,1,ppn);
            else
                svm = runSVMnoK(d,p,1,0);
            end
        else
%             svm = runSVMRepetitions4(d,p,1,ppn);
            svm = runSVMRepetitions5(d,p,1,ppn);
        end
        svm.r.dur = toc(stime);
        display(['This took ' datestr(1/24/60/60*svm.r.dur, 'HH:MM:SS') 'sec'])
        
        %         keyboard
        %
        %         length(find(svm.r.hitRate~=svm.r.hitRate_d))/numel(svm.r.hitRate)
        %         length(find(svm.r.hitRate~=svm.r.hitRate_pat))/numel(svm.r.hitRate)
        %         length(find(svm.r.hitRate~=svm.r.hitRate_patRel))/numel(svm.r.hitRate)
        %
        %
        %
        %         close all
        %         figure('WindowStyle','docked')
        %         errorbar(squeeze(mean(svm.r.pat,1)),squeeze(std(svm.r.pat,[],1)))
        %         legend({'pat1','pat2'}); ylim([-1 1])
        %         title('pat')
        %         figure('WindowStyle','docked')
        %         errorbar(squeeze(mean(svm.r.patRel,1)),squeeze(std(svm.r.patRel,[],1)))
        %         legend({'pat1','pat2'}); ylim([-1 1])
        %         title('patRel')
        %
        %         pat = squeeze(mean(svm.r.pat,1));
        %         figure('WindowStyle','docked')
        %         bar([mean(pat(1:end/2,:),1);mean(pat(end/2+1:end,:),1)]); hold on
        %         offset = 0.15;
        %         errorbar([1-offset 1+offset;2-offset 2+offset],[mean(pat(1:end/2,:),1);mean(pat(end/2+1:end,:),1)],[std(pat(1:end/2,:),[],1);std(pat(end/2+1:end,:),[],1)],'.')
        %         legend({'pat1','pat2'})
        %         title('pat')
        %
        %         pat = squeeze(mean(svm.r.patRel,1));
        %         figure('WindowStyle','docked')
        %         bar([mean(pat(1:end/2,:),1);mean(pat(end/2+1:end,:),1)]); hold on
        %         offset = 0.15;
        %         errorbar([1-offset 1+offset;2-offset 2+offset],[mean(pat(1:end/2,:),1);mean(pat(end/2+1:end,:),1)],[std(pat(1:end/2,:),[],1);std(pat(end/2+1:end,:),[],1)],'.')
        %         legend({'pat1','pat2'})
        %         title('patRel')
        
        
        %% Output and plot
        if stopThis
            break
        else
            if strcmp(p.k,'none')
                save(fullfile(dataDirOut,[dataFileOut '.mat']),'svm','d','p','perm','md','Fmask','-v7.3')
                continue
            end
            %reduce size and save
            svmOrig = svm;
            svm.w = [];
            svm.a = [];
            svm.amp = [];
            svm.delay = [];
            svm.xValCorr = [];
            save(fullfile(dataDirOut,[dataFileOut '.mat']),'svm','d','p','perm','md','Fmask')
            display(['Data saved as ' dataFileOut])
            
            %plot
            switch svm.p.algorithm
                case 'runSVM_RFE'
                    [h,~] = plotSVM(svm);
                    saveas(h,fullfile(dataDirOut,[dataFileOut '.jpg']),'jpg')
                    close(h)
                    
                    %                     [h,~] = plotSVM2(svm,[],[],[],[],'actualBinomDistThresh',plotIt);
                    %                     saveas(h,fullfile(dataDirOut,[dataFileOut '_thresh.jpg']),'jpg')
                    %                     close(h)
                case 'runSVM'
                    [h,~] = plotSVM2(svm,[],[],[],[],[],plotIt);
                    saveas(h,fullfile(dataDirOut,[dataFileOut '.jpg']),'jpg')
                    close(h)
                    
                    [h,~] = plotSVM2(svm,[],[],[],[],'actualBinomDistThresh',plotIt);
                    saveas(h,fullfile(dataDirOut,[dataFileOut '_thresh.jpg']),'jpg')
                    close(h)
                case 'runSVM_patAct2'
                    error('double-check that')
                    [h,~] = plotSVM_patAct(svm);
                    saveas(h,fullfile(dataDirOut,[dataFileOut '.jpg']),'jpg')
                    close(h)
                    
                case {'runSVM_patAct4','runSVM_patAct6'}
                    if svm.p.doCrossInfo
                        [h,~] = plotSVM_patAct(svm);
                    else
                        [h,~] = plotSVM_patAct5(svm);
                    end
                    saveas(h,fullfile(dataDirOut,[dataFileOut '.jpg']),'jpg')
                    close(h)
                case 'runSVMOneRepetition'
                    h = plotSVM_forSVMOneRepetition(svm);
                    saveas(h,fullfile(dataDirOut,[dataFileOut '.jpg']),'jpg')
                    close(h)
                case {'runSVM_patAct5'}
                    [h,~] = plotSVM_patAct5(svm);
                    saveas(h,fullfile(dataDirOut,[dataFileOut '.jpg']),'jpg')
                    close(h)
                case {'runSVM_patAct7'}
                    [h,~] = plotSVM_patAct7(svm);
                    saveas(h,fullfile(dataDirOut,[dataFileOut '.jpg']),'jpg')
                    close(h)
                otherwise
                    error('X')
            end
        end
    else
        
        %% Load already performed SVM analysis
        load(allFile{taskInd},'p','svm','d','perm')
        if ~exist('plotIt','var')
            plotIt = 0;
        end
        switch svm.p.algorithm
            case {'runSVM','runSVM_patAct6','runSVMOneRepetition'}
                if exist('altPerm','var')
                    if isfield(svm,'rP')
                        if strcmp(svm.pP.algorithm,'runSVMOneRepetition')
                            if size(svm.rP.newRes.hitRate,1)>=altPerm.numPerm
                                perm.doPerm = 0;
                            else
                                perm.doPerm = 1;
                            end
                        else
                            if size(svm.rP.hitRate,1)>=altPerm.numPerm;
                                perm.doPerm = 0;
                            else
                                perm.doPerm = 1;
                            end
                        end
                    else
                        perm.doPerm = 1;
                    end
                    plotIt = 1;
                    allFields = fields(altPerm);
                    for i = 1:length(allFields)
                        curField = allFields{i};
                        perm.(curField) = altPerm.(curField);
                    end
                end
            case {'runSVM_RFE','runSVM_patAct4'}
                %             case 'runSVMOneRepetition'
                
            otherwise
                error('did not specify what to do in that case')
        end
        
        
        %% Plot it
        if plotIt && ~strcmp(p.k,'none')
            switch svm.p.algorithm
                case 'runSVM'
                    [h1,h2] = plotSVM2(svm);
                case 'runSVM_RFE'
                    [h1,h2] = plotSVM(svm);
                case {'runSVM_patAct4','runSVM_patAct6'}
                    if svm.p.doCrossInfo
                        [h,~] = plotSVM_patAct(svm);
                    else
                        [h,~] = plotSVM_patAct5(svm);
                    end
                case 'runSVMOneRepetition'
                    h = plotSVM_forSVMOneRepetition(svm);
            end
        end
    end
    
    
    %% Permutation test
    if perm.doPerm
        %
        %         %% House keeping
        %         %Define feature levels
        %         if strcmp(p.algorithm,'SVM_RFE')
        %             error('code that')
        %             altFeatLev = p.featLevList;
        %         else
        %             if perm.doPermAtMax
        %                 altFeatLev =[];
        %                 [~,b] = max(mean(svm.r.hitRate,1));
        %                 altFeatLev =[altFeatLev p.featLevList(b)];
        %                 [~,b] = min(mean(svm.r.hitRate,1));
        %                 altFeatLev =[altFeatLev p.featLevList(b)];
        %             elseif perm.doPermAtAll
        %                 altFeatLev = p.featLevList;
        %             else
        %                 error('code that')
        %                 if p.nFeatLev>3
        %                     altFeatLev=[]; % define manually which point to run permutation test on
        %                     %                 altFeatLev=[9]; % define manually which point to run permutation test on
        %                     if stopIt
        %                         keyboard
        %                         %                 altFeatLev=[]; % define manually which point to run permutation test on
        %                     else
        %                         %                     keyboard
        %                     end
        %                 else
        %                     altFeatLev=p.featLevList;
        %                 end
        %             end
        %         end
        %
        %
        %Prepare p
        if isfield(svm,'pP')
            p = svm.pP;
            if exist('altPerm','var') && isfield(altPerm,'compileEachNperm')
                p.repeat = altPerm.compileEachNperm;
            end
            %             if exist('altFeatLev','var')
            %                 if ~isempty(altFeatLev)
            %                     p.altFeatLev = altFeatLev;
            %                 else
            %                     altFeatLev = p.altFeatLev;
            %                 end
            %             else
            %                 error('Must specify this because lazy programmer')
            %             end
        else
            svm.pP = svm.p;
            svm.pP.allPerm = nan(1,length(d.label));
            svm.pP.numPerm = perm.numPerm;
            svm.pP.repeat = perm.compileEachNperm;
            p = svm.pP;
            %             if exist('altFeatLev','var')
            %                 if ~isempty(altFeatLev)
            %                     p.altFeatLev = altFeatLev;
            %                 else
            %                     p.altFeatLev = p.featLevList;
            %                     altFeatLev = p.altFeatLev;
            %                 end
            %             else
            %                 error('Must specify this because lazy programmer')
            %             end
        end
        
        %Prepare r
        if isfield(svm,'rP')
            r = svm.rP;
        else
            if isfield(svm.r,'newRes')
                if strcmp(p.algorithm,'runSVMOneRepetition')
                    svm.rP.newRes = svm.r.newRes;
                    allFields = fields(svm.rP.newRes);
                    for i = 1:length(allFields)
                        switch class(svm.rP.newRes.(allFields{i}))
                            case 'double'
                                dim = size(svm.rP.newRes.(allFields{i})); dim(1) = 0;
                                svm.rP.newRes.(allFields{i}) = nan(dim);
                            case 'cell'
                                for ii = 1:size(svm.rP.newRes.(allFields{i}),2)
                                    switch class(svm.rP.newRes.(allFields{i}){1,ii})
                                        case 'double'
                                            dim = size(svm.rP.newRes.(allFields{i}){1,ii}); dim(1) = 0;
                                            svm.rP.newRes.(allFields{i}){1,ii} = nan(dim);
                                        case 'char'
                                            svm.rP.newRes.(allFields{i}){1,ii} = '';
                                        otherwise
                                            error('X')
                                    end
                                end
                                if length(allFields{i})>=4 && strcmp(allFields{i}(1:4),'info')
                                    svm.rP.newRes.(allFields{i})(2:end,:) = [];
                                end
                            otherwise
                                error('X')
                        end
                    end
                    
                else
                    svm.rP.hitRate = [];
                    svm.rP.d = [];
                    svm.rP.w = [];
                    svm.rP.a = [];
                    svm.rP.svmBias = [];
                end
            else
                svm.rP.hitRate = [];
                svm.rP.hitRateComb = [];
                svm.rP.d = [];
                svm.rP.w = [];
                svm.rP.a = [];
            end
            r = svm.rP;
        end
        
        
        
        %         % Compute only at all level of feature selection if not specified
        %         if ~exist('altFeatLev','var') || isempty(altFeatLev)
        %             altFeatLev = p.featLevList;
        %         end
        
        
        
        % Loop over permutations in the bloc
        i=0;
        permDone = false;
        allHitRate = [];
        while i<perm.numPermPerBloc && ~permDone
            % Run some SVMs with independent permutations (permutation done
            % within runSVM5.m within runSVMRepetitions4.m)
            if strcmp(p.k,'none')
                cursvm = runSVMnoK(d,p,0,ppn,1);
            else
                cursvm = runSVMRepetitions4(d,p,0,ppn,1);
            end
            if strcmp(cursvm.p.k,'none')
                curRes.w = cursvm.r.crossInfo.w;
                curRes.a = cursvm.r.crossInfo.a;
                curRes.info = cursvm.r.crossInfo_dim1Train;
                
                % Compile permutations in bloc
                i = i+size(curRes.w,1);
                
                % Compile permutations total
                %                     r.hitRate = cat(1,r.hitRate,curRes.hitRate);
                %                     r.d = cat(1,r.d,curRes.d);
                %                     r.svmBias = cat(1,r.svmBias,curRes.svmBias);
                r.w = cat(1,r.w,curRes.w);
                r.a = cat(1,r.a,curRes.a);
                r.info = curRes.info;
                if size(r.w,1)>=perm.numPerm
                    permDone = true;
                end
            else
                if isfield(cursvm.r,'newRes')
                    if strcmp(cursvm.p.algorithm,'runSVMOneRepetition')
                        
                        allFields = fields(cursvm.r.newRes);
                        for ii = 1:length(allFields)
                            switch class(cursvm.r.newRes.(allFields{ii}))
                                case 'double'
                                    r.newRes.(allFields{ii}) = cat(1,r.newRes.(allFields{ii}),cursvm.r.newRes.(allFields{ii}));
                                case 'cell'
                                    for iii = 1:size(cursvm.r.newRes.(allFields{ii}),2)
                                        switch class(cursvm.r.newRes.(allFields{ii}){1,iii})
                                            case 'double'
                                                r.newRes.(allFields{ii}){1,iii} = cat(1,r.newRes.(allFields{ii}){1,iii},cursvm.r.newRes.(allFields{ii}){1,iii});
                                            case 'char'
                                                r.newRes.(allFields{ii}){1,iii} = cursvm.r.newRes.(allFields{ii}){1,iii};
                                            otherwise
                                                error('X')
                                        end
                                    end
                                otherwise
                                    error('X')
                            end
                        end
                        
                        % Compile permutations in bloc
                        i = i+size(cursvm.r.newRes.hitRate,1);
                        if size(cursvm.r.newRes.hitRate,1)>=perm.numPerm
                            permDone = true;
                        end
                    else
                        
                        
                        curRes.hitRate = cursvm.r.newRes.hitRate;
                        curRes.d = cursvm.r.newRes.d;
                        curRes.svmBias = cursvm.r.newRes.svmBias;
                        curRes.w = squeeze(mean(cursvm.r.crossInfo.w,3));
                        curRes.a = squeeze(mean(cursvm.r.crossInfo.a,3));
                        curRes.info = cursvm.r.newRes.info(1,:);
                        
                        % Compile permutations in bloc
                        i = i+size(curRes.hitRate,1);
                        
                        % Compile permutations total
                        r.hitRate = cat(1,r.hitRate,curRes.hitRate);
                        r.d = cat(1,r.d,curRes.d);
                        r.svmBias = cat(1,r.svmBias,curRes.svmBias);
                        r.w = cat(1,r.w,curRes.w);
                        r.a = cat(1,r.a,curRes.a);
                        r.info = curRes.info;
                        if size(r.hitRate,1)>=perm.numPerm
                            permDone = true;
                        end
                    end
                else
                    cursvm.r.hitRate = nanmean(cursvm.r.hitRate,4);
                    cursvm.r.hitRateComb = cursvm.r.combInfo.hitRate;
                    tmp = nan(size(cursvm.r.crossInfo.d,1),size(cursvm.r.crossInfo.d,2),size(cursvm.r.crossInfo.d,4));
                    for ii = 1:size(cursvm.r.crossInfo.d,1)
                        for iii = 1:size(cursvm.r.crossInfo.d,4)
                            tmp(ii,:,iii) = diag(squeeze(cursvm.r.crossInfo.d(ii,:,:,iii)));
                        end
                    end
                    
                    cursvm.r.crossInfo.d = tmp;
                    cursvm.r.crossInfo.w = nanmean(cursvm.r.crossInfo.w,3);
                    cursvm.r.crossInfo.a = nanmean(cursvm.r.crossInfo.a,3);
                    
                    % Compile permutations in bloc
                    i = i+size(cursvm.r.hitRate,1);
                    
                    % Compile permutations total
                    r.hitRate = cat(1,r.hitRate,cursvm.r.hitRate);
                    r.hitRateComb = cat(1,r.hitRateComb,cursvm.r.hitRateComb');
                    r.d = cat(1,r.d,cursvm.r.crossInfo.d);
                    r.w = cat(1,r.w,squeeze(cursvm.r.crossInfo.w));
                    r.a = cat(1,r.a,squeeze(cursvm.r.crossInfo.a));
                    if size(r.hitRate,1)>=perm.numPerm
                        permDone = true;
                    end
                end
            end
            % Put data back in the main structure
            svm.pP = p;
            svm.rP = r;
            
            % Compute non-param thresh so far
            if ~strcmp(p.k,'none')
                if ~strcmp(svm.pP.algorithm,'runSVMOneRepetition')
                    svm.rP.thresh = prctile(svm.rP.hitRate,[2.5 5 95 97.5],1);
                end
            end
            %             if strcmp(svm.pP.infoComb,'layeredSVM')
            %                 svm.rP.threshComb = prctile(svm.rP.hitRateComb,[2.5 5 95 97.5],1);
            %             end
            
            % Plot
            switch p.algorithm
                case {'runSVM_patAct6'}
                    if ~strcmp(svm.pP.k,'none')
                        svm.plotIt = plotIt;
                        if ~exist('h','var')
                            h = [];
                        end
                        if svm.p.doCrossInfo
                            error('double-check that')
                            [h,~] = plotSVM_patAct(svm);
                        else
                            [h,~] = plotSVM_patAct5(svm,h,1,1);
                            saveas(h,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '.jpg']),'jpg')
                            %                         close(h)
                        end
                    end
                case 'runSVMOneRepetition'
                    h = plotSVM_forSVMOneRepetition(svm,h,1);
                    saveas(h,fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '.jpg']),'jpg')
                    
                otherwise
                    error('double-check that')
                    plotNameList = {'actualMean' 'negThresh' 'posThresh' 'permMean' 'actualThresh'};
                    svm.plotIt = plotIt;
                    if ~exist('h','var')
                        h = zeros(1,length(plotNameList));
                        if exist('h1','var')
                            h(1) = h1;
                        end
                    end
                    for iii = 1:length(plotNameList)
                        if h(iii)
                            h(iii) = plotSVM2(svm,h(iii),[],[],[],plotNameList{iii});
                        else
                            h(iii) = plotSVM2(svm,[],[],[],[],plotNameList{iii});
                        end
                        saveas(h(iii),fullfile(p.dataDirOut,[p.dataFileOut '_' plotNameList{iii} '.jpg']),'jpg')
                    end
                    
            end
            drawnow
            
            
            % Save
            if ~exist(fullfile(p.dataDirOut,[p.dataFileOut '.mat']),'file')
                p.dataDirOut = altMainDir;
            end
            save(fullfile(p.dataDirOut,[p.dataFileOut '.mat']),'svm','-append','-v7.3')
            
            
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
            %             if svm.p.doSVM
            %                 allHitRate = nan(perm.compileEachNperm,length(altFeatLev),length(altFeatLev));
            %             end
            % %             if svm.p.doCrossVal
            % %                 allRsquared = nan(perm.compileEachNperm,length(altFeatLev),length(altFeatLev));
            % %             end
            %
            %             for ii = 1:perm.compileEachNperm
            %                 %% Check if we have enough
            %                 if svm.p.doSVM
            %                     permInBloc = findNperm(allHitRate);
            %                     permTotal = findNperm(r.hitRate);
            %                 else
            %                     permInBloc = findNperm(allRsquared);
            %                     permTotal = findNperm(xValCorr.Rsquared);
            %                 end
            %                 nPerm = permTotal + permInBloc;
            %                 if all(nPerm>=perm.numPerm)
            %                     permDone = 1;
            %                     break
            %                 end
            %
            %                 %% Do it
            %                 display(['Perm ' num2str(i+1) '/' num2str(perm.numPermPerBloc) ' in bloc']);
            %                 display(['Perm ' num2str(max(permTotal)+1) '/' num2str(perm.numPerm) ' total']);
            %                 i = i+1;
            %
            %                 tic
            %                 curd = d; tmpp = p;
            %
            %                 % Shuffle data
            %                 %randomly swap pairs
            % %                 for obs = 1:size(curd.xData,1)/2
            % %                     ind = [obs obs+size(curd.xData,1)/2]; indShuf = ind;
            % %                     if round(rand(1,1)); indShuf = fliplr(indShuf); end
            % %                     curd.xData(ind,:) = d.xData(indShuf,:);
            % %                 end
            %                 %shuffle the whole thing
            %                 curd.xData = d.xData(randperm(size(d.xData,1)),:);
            %
            %
            %                 tmpp.repeat = 1; %only one repeat for permutation test
            %                 tmpp.featLevList = altFeatLev;
            %                 tmpp.nFeatLev = length(tmpp.featLevList);
            %                 if ~all(ismember(tmpp.featLevList,svm.p.featLevList))
            %                     error('Must feature levels for permutation that were used for actual test')
            %                 end
            %
            %                 cursvm = runSVMRepetitions3(curd,tmpp,0,0);
            %                 allHitRate(ii,:,:) = mean(cursvm.r.hitRate,1);
            %
            %                 dur(ii) = toc;
            %                 display(['Took ' num2str(round(dur(ii))) 'sec']);
            %
            %             end
            %
            %             %% Compile
            %
            %             r.hitRate = cat(1,r.hitRate,allHitRate);
            % %             r.hitRate = compilePermRes(r.hitRate,svm.p.featLevList,svm.p.featLevList,allHitRate,altFeatLev,altFeatLev);
            % %             xValCorr.Rsquared = compilePermRes(xValCorr.Rsquared,svm.p.featLevList,allRsquared,altFeatLev);
            %             if ~isfield(r,'dur')
            %                 r.dur = dur';
            %             else
            %                 if ~exist('dur','var')
            %                     dur=0;
            %                 end
            %                 r.dur = [r.dur; dur'];
            %             end
            %
            %             % Put data back in the main structure
            %             svm.pP = p;
            %             svm.rP = r;
            %
            %             % Compute non-param thresh so far
            %             svm.rP.thresh = prctile(svm.rP.hitRate,[2.5 5 95 97.5],1);
            %
            %             % Plot
            %             plotNameList = {'actualMean' 'negThresh' 'posThresh' 'actualThresh' 'permMean'};
            %             svm.plotIt = plotIt;
            %             if ~exist('h','var')
            %                 h = zeros(1,length(plotNameList));
            %                 if exist('h1','var')
            %                     h(1) = h1;
            %                 end
            %             end
            %             for iii = 1:length(plotNameList)
            %                 if h(iii)
            %                     h(iii) = plotSVM2(svm,h(iii),[],[],[],plotNameList{iii});
            %                 else
            %                     h(iii) = plotSVM2(svm,[],[],[],[],plotNameList{iii});
            %                 end
            %                 saveas(h(iii),fullfile(p.dataDirOut,[p.dataFileOut '_' plotNameList{iii} '.jpg']),'jpg')
            %             end
            %
            %
            %             % Save
            %             if ~exist(fullfile(p.dataDirOut,[p.dataFileOut '.mat']),'file')
            %                 p.dataDirOut = altMainDir;
            %             end
            %             save(fullfile(p.dataDirOut,[p.dataFileOut '.mat']),'svm','-append')
        end
    end
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
end
