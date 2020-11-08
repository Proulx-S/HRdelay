function [ taskDone ] = runTask8(allFile,perm,altMainDir,guillimin,ppn)

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
    clearvars -except taskInd allFile plotIt doPerm dataDirOut perm altMainDir guillimin ppn taskDone
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
        
        p.nTrialPerRun = str2double(p.split(1));
        if p.splitHalf
            d.xData = d.xData(p.splitHalf:2:end,:);
            d.normData = d.normData(p.splitHalf:2:end,:);
            d.label = d.label(p.splitHalf:2:end,:);
            d.sessionLabel = d.sessionLabel(p.splitHalf:2:end,:);
            p.nObs = p.nObs/2;
            p.nTrialPerRun = p.nTrialPerRun/2;
        end
        d.runLabel = reshape(repmat(1:p.nObs /p.nTrialPerRun,p.nTrialPerRun,1),p.nObs,1);
        
        
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
        end
        svm.dur = toc(stime);
        display(['This took ' datestr(1/24/60/60*svm.dur, 'HH:MM:SS') 'sec'])
        
        
        
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
            
            %Plot
%             h = plotSVM_forParaDis(svm);
%             h = plotSVM_forSVMOneRepetition2(svm);
%             h = plotSVM_forSVMOneRepetition3(svm);
            h = plotSVM_forParaAcc(svm,[],[],0);
            for anaInd = 1:length(h)
                saveas(h(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '.jpg']),'jpg')
                close(h(anaInd))
            end
            h2 = plotSVM_forParaAcc(svm,[],[],1);
            for anaInd = 1:length(h2)
                saveas(h2(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_abs.jpg']),'jpg')
                close(h2(anaInd))
            end
            h3 = plotSVM_forParaDist(svm,[],[],1);
            for anaInd = 1:length(h3)
                saveas(h3(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_distAbs.jpg']),'jpg')
                close(h3(anaInd))
            end
        end
    else
        %% Load already performed SVM analysis
        load(allFile{taskInd},'svm','md','Fmask')
        h = plotSVM_forParaAcc(svm,[],[],0);
        h2 = plotSVM_forParaAcc(svm,[],[],1);
        h3 = plotSVM_forParaDist(svm,[],[],1);
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
            if isfield(svm.perm,'rTe')
                svm.perm.rTe = clearAllDataField2(svm.perm.rTe,[],1);
            end
        end
%         svmCompiled = svm;
%         svm = svmCompiled.perm;
        
        
        % Loop over permutations in the bloc
        if length(svm.perm.rTe{1}.acc)<p.perm.numPerm
            permDone = false;
        else
            permDone = true;
%             %Plot
%             h = plotSVM_forSVMOneRepetition3(svm,h); % for all in different plots    
        end
        
        i=0;
        while i<svm.p.perm.numPermPerBloc && ~permDone
            svm_cur = runSVMRepetitions5(d,p,1,ppn);

            tic
            display(['Compiling and saving as ' svm.p.dataFileOut])
            % Compile permutations in bloc
            svm_cur.rTr = clearAllDataField(svm_cur.rTr,[],0);
            svm_cur.rTe = clearAllDataField(svm_cur.rTe,[],0);
            
            svm.perm.rTr = compileRep(svm.perm.rTr,svm_cur.rTr);
            svm.perm.rTe = compileRep(svm.perm.rTe,svm_cur.rTe);
            
            i = i+length(svm_cur.rTe{1}.acc);
            if length(svm.perm.rTe{1}.acc)>=p.perm.numPerm
                permDone = true;
            end
            
            %Plot
%             h = plotSVM_forParaAcc(svm,h);
%             h = plotSVM_forParaDis(svm,h);
%             h = plotSVM_forSVMOneRepetition2(svm,h);
%             h = plotSVM_forSVMOneRepetition3(svm,h); % for all in different plots
            h = plotSVM_forParaAcc(svm,h,[],0);
            for anaInd = 1:length(h)
                saveas(h(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '.jpg']),'jpg')
            end
            h2 = plotSVM_forParaAcc(svm,h2,[],0);
            for anaInd = 1:length(h2)
                saveas(h2(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_abs.jpg']),'jpg')
            end
            h3 = plotSVM_forParaAcc(svm,h3,[],0);
            for anaInd = 1:length(h3)
                saveas(h3(anaInd),fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '__' svm.p.infoComb{anaInd} '_distAbs.jpg']),'jpg')
            end
            
            %Save
            try
                save(fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '.mat']),'svm','-append')
            catch
                save(fullfile(svm.p.dataDirOut,[svm.p.dataFileOut '.mat']),'svm','-append','-v7.3')
            end
            display(['Done in'])
            toc
            
        end
        for anaInd = 1:length(h)
            close(h(anaInd))
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
