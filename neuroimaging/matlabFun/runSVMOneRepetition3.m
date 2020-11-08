function [rTr,rVal,rTe] = runSVMOneRepetition3(d,p,verbose,rep,tmpStrResp)
if ~exist('tmpStrResp','var')
    global tmpStrResp
end

if ~exist('verbose','var')
    verbose = 0;
end
% if ~exist('doPerm','var')
%     doPerm = 0;
% end

%% K-folding and permutations
[d.crossVal,d.label] = defineKnPerm(d,p);


%% Multilevel svm
if ~iscell(p.infoComb)
    p.infoComb = {p.infoComb};
end

dOrig = d;
for infoCombInd = 1:length(p.infoComb)
    d = dOrig;
%     tmpStrResp = {'respSVM__pcaNone_pcaNone','respSVM__pcaRun_pcaNone','respSVM__pcaTime_pcaNone','respSVM__pcaNone_pcaSpace','respSVM__pcaRun_pcaSpace','respSVM__pcaTime_pcaSpace','respSVM__pcaNone_pcaRun','respSVM__pcaRun_pcaRun','respSVM__pcaTime_pcaRun'};
    if ismember(p.infoComb(infoCombInd),tmpStrResp)
        d.xData = d.respData;
        d.normData = d.respDataNorm;
    else
        if isfield(d,'respData')
            d = rmfield(d,'respData');
        end
    end
    curP = p; curP.infoComb = p.infoComb{infoCombInd};
    if ~isnumeric(p.cList)
        curP.C.subj = curP.C.subj(strmatch(curP.subj,curP.C.subj));
        curP.C.Cmax = curP.C.Cmax(strmatch(curP.subj,curP.C.subj),:);
    else
        if length(p.cList)==length(p.infoComb)
            curP.C = p.cList(infoCombInd);
        elseif length(p.cList)==1
            curP.C = p.cList;
        end
    end
    switch curP.infoComb
%         case {'all','all_cross1','all_cross2','ALrev','ALrev_cross1','ALrev_cross2'}
%             tmpStr = strsplit(curP.infoComb,'_'); if strcmp(tmpStr{end}(1:5),'cross'); curP.subCrossVal = str2num(tmpStr{end}(6)); else curP.subCrossVal = 0; end
%             [rTr{infoCombInd},rTe{infoCombInd}] = runMultiLevelSVM2(d,curP);

%         case {'subFold_ALrev','subFold_noALrev','voxComb_subFold','infoComb_subFold','catComb_subFold'}
        case {'lowLevel_subFold','voxComb_subFold','infoComb_subFold','catComb_subFold' 'lowLevel_subFold_AL','voxComb_subFold_AL','infoComb_subFold_AL','catComb_subFold_AL','infoComb_CsubFold'}
            [rTr{infoCombInd},rVal{infoCombInd},rTe{infoCombInd}] = runMultiLevelSVM6(d,curP,1);
            
%             delayDiff = wrapToPi(angle(d.xData(d.label==1,:))-angle(d.xData(d.label==2,:)));
%             delayDiff = circ_mean(delayDiff,[],1);
%             ampDiff = abs(d.xData(d.label==1,:))-abs(d.xData(d.label==2,:));
%             ampDiff = mean(ampDiff,1);
%             
%             hist(ampDiff(ampDiff>1 | ampDiff<-1))
%             hist(ampDiff)
%             
%             scatter(delayDiff,ampDiff)
%             hist(angle(d.xData(1,:)))
%             [rho, pval] = circ_corrcl(delayDiff,ampDiff);
            
        case {'lowLevel','voxComb','infoComb','catComb','infoComb_Copt','voxComb_Copt'}
            [rTr{infoCombInd},rTe{infoCombInd}] = runMultiLevelSVM7(d,curP,1);
            rVal = [];
        case {'pcaOnResp_lowLevel','pcaOnTS_lowLevel'}
%             pcaRes = [];
            coeff = [];
            score = [];
            latent = [];
            xData = d.xData; % run x vox x time
            for vox = 1:size(xData,2)
                curData = permute(xData(:,vox,:),[3 1 2]); % time x run
                [coeff(:,:,vox),score(:,:,vox),latent(:,vox),~,~,~] = pca(curData);
            end
            d.xData = permute(coeff,[1 3 2]); % run x vox x pc <-- run x pc x vox
            [rTr{infoCombInd},rTe{infoCombInd}] = runPCASVM(d,curP,1);
            rVal = [];
%             rTr{infoCombInd}.pca = pcaRes;
        case 'respSVM__pcaNone_pcaNone'
            [rTr{infoCombInd},rTe{infoCombInd}] = runPCASVM(d,curP,1);
            rVal = [];           
        case 'respSVM__pcaRun_pcaNone'
            pcaRes = [];
            coeff = [];
            score = [];
            latent = [];
            tsquared = [];
            mu = [];
            xData = d.xData; % run x vox x time
            for vox = 1:size(xData,2)
                curData = permute(xData(:,vox,:),[3 1 2]); % time x run
                [coeff(:,:,vox),score(:,:,vox),latent(:,vox),tsquared(:,vox),~,mu(:,vox)] = pca(curData,'Centered',true); % coeff: run x pc x vox
            end
            d.xData = permute(coeff,[1 3 2]); % run x vox x pc <-- run x pc x vox
            
            [rTr{infoCombInd},rTe{infoCombInd}] = runPCASVM(d,curP,1);
            rVal = [];
            
            pcaRes.coeff = coeff;
            pcaRes.score = score;
            pcaRes.latent = latent;
            pcaRes.tsquared = tsquared;
            pcaRes.mu = mu;
            rTr{infoCombInd}.pca = pcaRes;

        case 'respSVM__pcaTime_pcaNone'
            pcaRes = [];
            coeff = [];
            score = [];
            latent = [];
            tsquared = [];
            mu = [];
            xData = d.xData; % run x vox x time
            for vox = 1:size(xData,2)
                curData = squeeze(xData(:,vox,:)); % run x time
                [coeff(:,:,vox),score(:,:,vox),latent(:,vox),tsquared(:,vox),~,mu(:,vox)] = pca(curData,'Centered',true); % coeff: time x pc x vox; score: run x pc x vox
            end
            d.xData = permute(score,[1 3 2]); % run x vox x pc
            
            [rTr{infoCombInd},rTe{infoCombInd}] = runPCASVM(d,curP,1);
            rVal = [];
            
            pcaRes.coeff = coeff;
            pcaRes.score = score;
            pcaRes.latent = latent;
            pcaRes.tsquared = tsquared;
            pcaRes.mu = mu;
            rTr{infoCombInd}.pca = pcaRes;

        case 'respSVM__pcaNone_pcaSpace'
            coeff = [];
            score = [];
            latent = [];
            mu = [];
            xData = d.xData; % run x vox x time
            for timeInd = 1:size(xData,3)
                curData = xData(:,:,timeInd); % run x vox
                [coeff(:,:,timeInd),score(:,:,timeInd),latent(:,timeInd),~,~,mu(:,vox)] = pca(curData); % coeff: vox x pc x time; score: run x pc x time
            end
            d.xData = score; % run x pc x time
            
            [rTr{infoCombInd},rTe{infoCombInd}] = runPCASVM(d,curP,1);
            rVal = [];
        case 'respSVM__pcaRun_pcaSpace'
            coeff = [];
            score = [];
            latent = [];
            mu = [];
            xData = d.xData; % run x vox x time
            for vox = 1:size(xData,2)
                curData = permute(xData(:,vox,:),[3 1 2]); % time x run
                [coeff(:,:,vox),score(:,:,vox),latent(:,vox),~,~,mu(:,vox)] = pca(curData); % coeff: run x pc x vox
            end
            d.xData = permute(coeff,[1 3 2]); % run x vox x pc1
            
            coeff = [];
            score = [];
            latent = [];
            mu = [];
            xData = d.xData; % run x vox x pcL1
            for pcInd = 1:size(xData,3)
                curData = xData(:,:,pcInd); % run x vox
                [coeff(:,:,pcInd),score(:,:,pcInd),latent(:,pcInd),~,~,mu(:,vox)] = pca(curData); % coeff: vox x pcL2 x pcL1; score: run x pcL2 x pcL1
            end
            d.xData = score; % run x pcL2 x pcL1
            
            [rTr{infoCombInd},rTe{infoCombInd}] = runPCASVM(d,curP,1);
            rVal = [];
            
        case 'respSVM__pcaTime_pcaSpace'
            coeff = [];
            score = [];
            latent = [];
            mu = [];
            xData = d.xData; % run x vox x time
            for vox = 1:size(xData,2)
                curData = squeeze(xData(:,vox,:)); % run x time
                [coeff(:,:,vox),score(:,:,vox),latent(:,vox),~,~,mu(:,vox)] = pca(curData); % coeff: time x pc x vox; score: run x pc x vox
            end
            d.xData = permute(score,[1 3 2]); % run x vox x pcL1
            
            coeff = [];
            score = [];
            latent = [];
            mu = [];
            xData = d.xData; % run x vox x pcL1
            for pcInd = 1:size(xData,3)
                curData = xData(:,:,pcInd); % run x vox
                [coeff(:,:,pcInd),score(:,:,pcInd),latent(:,pcInd),~,~,mu(:,vox)] = pca(curData); % coeff: vox x pcL2 x pcL1; score: run x pcL2 x pcL1
            end
            d.xData = score; % run x pcL2 x pcL1
            
            [rTr{infoCombInd},rTe{infoCombInd}] = runPCASVM(d,curP,1);
            rVal = [];
        case 'respSVM__pcaNone_pcaRun'
            coeff = [];
            score = [];
            latent = [];
            mu = [];
            xData = permute(d.xData,[2 1 3]); % vox x run x time
            for timeInd = 1:size(xData,3)
                curData = xData(:,:,timeInd); % vox x run
                [coeff(:,:,timeInd),score(:,:,timeInd),latent(:,timeInd),~,~,mu(:,vox)] = pca(curData); % coeff: run x pc x time; score: vox x pc x time
            end
            d.xData = coeff; % run x pc x time
            
            [rTr{infoCombInd},rTe{infoCombInd}] = runPCASVM(d,curP,1);
            rVal = [];
        case 'respSVM__pcaRun_pcaRun'
            coeff = [];
            score = [];
            latent = [];
            mu = [];
            xData = d.xData; % run x vox x time
            for vox = 1:size(xData,2)
                curData = permute(xData(:,vox,:),[3 1 2]); % time x run
                [coeff(:,:,vox),score(:,:,vox),latent(:,vox),~,~,mu(:,vox)] = pca(curData); % coeff: run x pc x vox
            end
            d.xData = permute(coeff,[1 3 2]); % run x vox x pcL1
            
            coeff = [];
            score = [];
            latent = [];
            mu = [];
            xData = permute(d.xData,[2 1 3]); % vox x run x pcL1
            for pcInd = 1:size(xData,3)
                curData = xData(:,:,pcInd); % vox x run
                [coeff(:,:,pcInd),score(:,:,pcInd),latent(:,pcInd),~,~,mu(:,vox)] = pca(curData); % coeff: run x pcL2 x pcL1; score: vox x pcL2 x pcL1
            end
            d.xData = coeff; % run x pcL2 x pcL1
            
            [rTr{infoCombInd},rTe{infoCombInd}] = runPCASVM(d,curP,1);
            rVal = [];
        case 'respSVM__pcaTime_pcaRun'
            coeff = [];
            score = [];
            latent = [];
            mu = [];
            xData = d.xData; % run x vox x time
            for vox = 1:size(xData,2)
                curData = squeeze(xData(:,vox,:)); % run x time
                [coeff(:,:,vox),score(:,:,vox),latent(:,vox),~,~,mu(:,vox)] = pca(curData); % coeff: time x pc x vox; score: run x pc x vox
            end
            d.xData = permute(score,[1 3 2]); % run x vox x pcL1
            
            coeff = [];
            score = [];
            latent = [];
            mu = [];
            xData = permute(d.xData,[2 1 3]); % vox x run x pcL1
            for pcInd = 1:size(xData,3)
                curData = xData(:,:,pcInd); % vox x run
                [coeff(:,:,pcInd),score(:,:,pcInd),latent(:,pcInd),~,~,mu(:,vox)] = pca(curData); % coeff: run x pcL2 x pcL1; score: vox x pcL2 x pcL1
            end
            d.xData = coeff; % run x pcL2 x pcL1
            
            [rTr{infoCombInd},rTe{infoCombInd}] = runPCASVM(d,curP,1);
            rVal = [];
            
        case 'pcaOnResp_pcaOnSpace'
            keyboard
            pcaRes = [];
            xData = d.xData;
            for vox = 1:size(xData,2)
                curData = permute(xData(:,vox,:),[3 1 2]);
                [coeff(:,:,vox),score(:,:,vox),latent(:,vox),~,~,mu(:,vox)] = pca(curData); %coeff: run x pcTime x vox
            end
            for pcTimeInd = 1:size(coeff,2)
                curData = squeeze(coeff(:,pcTimeInd,:)); %run x vox
                [coeff2(:,:,pcTimeInd),score2(:,:,pcTimeInd),latent2(:,pcTimeInd),~,~,mu2(:,vox)] = pca(curData); %score2: run x pcSpace x pcTime
            end
            clear coeff score latent coeff2 latent2
            d.xData = score2; % run x pcSpace x pcTime
            [rTr{infoCombInd},rTe{infoCombInd}] = runPCASVM(d,curP,1);
            rVal = [];
%             rTr{infoCombInd}.pca = pcaRes;
            
        case 'crossValDelay'
%             curP.nBin = 20;
%             curP.histWidth = 2*pi;
            [N,crossValDelay,binCent,crossValDelay2,delay2,binCent2] = runCrossValDelay2(d,curP,0);
            rTr{infoCombInd}.diff.N = N;
            rTr{infoCombInd}.diff.binCent = binCent;
            rTr{infoCombInd}.diff.binCent2 = binCent2;
            rTr{infoCombInd}.diff.delay2 = delay2;
            rTe{infoCombInd}.diff.crossValDelay = crossValDelay;
            rTe{infoCombInd}.diff.binCent = binCent;
            rTe{infoCombInd}.diff.crossValDelay2 = crossValDelay2;
            rTe{infoCombInd}.diff.binCent2 = binCent2;
            
%             curP.nBin = 60;
%             curP.histWidth = 20;
%             [N,crossValDelay,binCent] = runCrossValDelay2(d,curP,0);
%             rTr{infoCombInd}.t.N = N;
%             rTr{infoCombInd}.t.binCent = binCent;
%             rTe{infoCombInd}.t.crossValDelay = crossValDelay;
%             rTe{infoCombInd}.t.binCent = binCent;
            
            rVal = [];
            
            
        case 'crossValDelayQ'
            [delay,crossValDelay,N,edges,binCent] = runCrossValDelayQ(d,curP,0);
            rTr{infoCombInd}.diff.binCent = binCent;
            rTr{infoCombInd}.diff.edges = edges;
            rTr{infoCombInd}.diff.N = N;
            rTr{infoCombInd}.diff.delay = delay;
            
            rTe{infoCombInd}.diff.binCent = binCent;
            rTe{infoCombInd}.diff.edges = edges;
            rTe{infoCombInd}.diff.delay = crossValDelay;
            
            rVal = [];
        case 'respPCA'
            %Do some classification
            [rTr{infoCombInd},rTe{infoCombInd}] = runRespPCA(d,curP,0);
%             plotPCA(rTr{infoCombInd},p,2,0)
%             plotPCA(rTr{infoCombInd},p,0,1)
%             figure('WindowStyle','docked');
%             rTr{infoCombInd}.latent
%             [a,b] = sort(rTe{infoCombInd}.distT,'descend')
%             
%             bar(sort(rTe{infoCombInd}.distT,'descend'))
%             set(gca,'xTick',(1:24));
%             set(gca,'xTickLabel',cellstr(num2str(b')));
%             b
            rVal = [];
        case 'respPCA_PCAalign'
            
            %Do some classification
            runRespPCA_PCAalign2(d,curP,0);
            rTr = {};
            rVal = {};
            rTe = {};
            runRespPCA_PCAalign(d,curP,0);
        case 'none'
            error('X')
        otherwise
            error('X')
    end
end



function [crossVal,label] = defineKnPerm(d,p)

switch p.KnPerm
    case 'Prun_KtrialWLabel'
        %Permutations
        if p.doPerm
            runLabel = unique(d.runLabel);
            runLabel = reshape(repmat(runLabel(randperm(length(runLabel))),1,p.nTrialPerRun)',p.nObs,1);
            [~,b] = sort(runLabel);
            label = d.label(b);
            clear runLabel b
        else
            label = d.label;
        end
        
        %K-folding
        crossVal = nan(length(d.label),1);
        crossVal_tmp = (1:length(d.label)/2)';
        crossVal(label==1) = crossVal_tmp(randperm(length(crossVal_tmp)));
        crossVal(label==2) = crossVal_tmp(randperm(length(crossVal_tmp)));
    case 'Ptrial_KtrialWLabel'
        %Permutations
        if p.doPerm
            label = d.label(randperm(length(d.label)));
        else
            label = d.label;
        end
        
        %K-folding
        crossVal = nan(length(d.label),1);
        crossVal_tmp = (1:length(d.label)/2)';
        crossVal(label==1) = crossVal_tmp(randperm(length(crossVal_tmp)));
        crossVal(label==2) = crossVal_tmp(randperm(length(crossVal_tmp)));
    case 'Prun_KrunWLabel'
        %Permutations
        if p.doPerm
            runLabel = unique(d.runLabel);
            runLabel = reshape(repmat(runLabel(randperm(length(runLabel))),1,p.nTrialPerRun)',p.nObs,1);
            [~,b] = sort(runLabel);
            label = d.label(b);
            clear runLabel b
        else
            label = d.label;
        end
        
        %K-folding
        runLabel = unique(d.runLabel);
        crossVal = nan(length(d.label),1);
        
        crossVal_tmp = (1:length(runLabel)/2)';
        crossVal_tmp = crossVal_tmp(randperm(length(crossVal_tmp)));
        crossVal_tmp = reshape(repmat(crossVal_tmp,1,p.nTrialPerRun)',p.nObs/2,1);
        crossVal(label==1) = crossVal_tmp;
        
        crossVal_tmp = (1:length(runLabel)/2)';
        crossVal_tmp = crossVal_tmp(randperm(length(crossVal_tmp)));
        crossVal_tmp = reshape(repmat(crossVal_tmp,1,p.nTrialPerRun)',p.nObs/2,1);
        crossVal(label==2) = crossVal_tmp;
    otherwise
        error('X')
end

