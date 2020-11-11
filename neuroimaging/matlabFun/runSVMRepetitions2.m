function [svm,xValCorr] = runSVMRepetitions2(d,p,verbose)
global err
if ~exist('verbose','var')
    verbose = 1;
end

% keyboard
% phase = angle(d.xData);
% phase_shift = circ_mean(phase,[],1);
% phase = phase-repmat(phase_shift,size(phase,1),1)+45/180*pi;
% [x,y] = pol2cart(phase,1);
% d.xData = complex(x,y);

if ischar(p.k)
    switch p.k
        case 'auto'
            k = p.nObs/2;
        case 'autoRun'
            if strcmp(p.split(2:end),'perRun')
                k = p.nObs/2/str2double(p.split(1)); % to keep k constant (at the equivalent of leave-one-run-pair-out) across different spliting of the data
            else
                error('')
            end
        otherwise
            error('')
    end
else
    k = p.k;
end
% if mod(length(d.label)/k,1)
%     error('number of data point is not a multiple of k (programmer in a rush)')
% end
svm.w = nan(p.nFeatLev,p.nFeatures,k,p.repeat);
svm.a = nan(p.nFeatLev,p.nFeatures,k,p.repeat);
svm.r.hitRate = nan(p.repeat,p.nFeatLev,p.nFeatLev);
if p.thresholdData; svm.mi = nan(k,p.nFeatures,p.repeat); end
svm.amp = nan(p.nFeatLev,p.nFeatures,k,p.repeat);
svm.delay = nan(p.nFeatLev,p.nFeatures,k,p.repeat);
if strcmp(p.keepInfo,'p_then_m') || strcmp(p.keepInfo,'m_then_p')
    svm.alt.w = nan(p.nFeatLev,p.nFeatures,k,p.repeat);
    svm.alt.a = nan(p.nFeatLev,p.nFeatures,k,p.repeat);
    svm.alt.r.hitRate = nan(p.repeat,p.nFeatLev);
    svm.alt.amp = nan(p.nFeatLev,p.nFeatures,k,p.repeat);
    svm.alt.delay = nan(p.nFeatLev,p.nFeatures,k,p.repeat);
end
xValCorr.Rsquared = nan(p.repeat,p.nFeatLev);
xValCorr.RHO = nan(p.repeat,p.nFeatLev);

% trAll = zscore(abs(d.xData(:,p.funcROI.vec.sortInd)),[],1);
% tmpSort = mean(trAll(1:size(trAll,1)/2,:),1)<mean(trAll(size(trAll,1)/2+1:end,:),1);
% trTmp = trAll;
% trTmp(1:size(trAll,1)/2,tmpSort) = trAll(size(trAll,1)/2+1:end,tmpSort);
% trTmp(size(trAll,1)/2+1:end,tmpSort) = trAll(1:size(trAll,1)/2,tmpSort);
% figure('windowStyle','docked');
% imagesc(trTmp,[-2 2])
% clear trAll tmpSort trTmp

% keyboard
% 
% clear tmp
% tmp(:,:,1) = abs(d.xData(1:end/2,:));
% tmp(:,:,2) = abs(d.xData(end/2+1:end,:));
% tmp(:,:,3) = abs(d.normData);
% tmp = reshape(tmp,[size(tmp,1)*size(tmp,3) size(tmp,2)])';
% tmp = corr(tmp);
% close all
% figure('windowStyle','docked');
% imagesc(tmp,[0 1])
% 
% clear tmp1 tmp2
% tmp1(:,:,1) = [1 2 3; 4 5 6];
% tmp1(:,:,2) = [7 8 9;10 11 12];
% tmp1 = reshape(tmp1,[size(tmp1,1)*size(tmp1,3) size(tmp1,2)]);
% tmp2 = cat(1,[1 2 3; 4 5 6],[7 8 9;10 11 12]);
% 
% tmp1 = cat(1,abs(d.xData),abs(d.normData));
% tmp1 = corr(tmp1');
% close all
% figure('windowStyle','docked');
% imagesc(tmp1)
% 
% 
% imagesc(corr(tmp,tmp1'))


% if p.regSession
%     % Regress out session effect voxel-wise, incorporating plaid
%     %     metric.ampRect =  zscore(metricOrig.ampRect,[],1);
%     d = regressOutSessionVoxWise(d,p,0);
%     %%%%%%%
%     % Might have to insert that later to remove testing data from this
%     % normalization process (but for now, at least permutation is done before)
% end

if isfield(d,'crossVal')
    genCrossVal = 0;
else
    genCrossVal = 1;
end

for rep = 1:p.repeat
    if verbose
        tic
        display(['Starting repetition ' num2str(rep) '/' num2str(p.repeat)])
    end
    
    
    
    
    % Define cross validation folding
    if genCrossVal
        if ischar(p.k)
            switch p.k
                case 'auto'
                    %Leave-one-pair-out, cond label balanced only:
                    %Random selection of two examples for the testing set at each fold.
                    %Each example goes into the testing set once and only once. Does not
                    %care about sessions.
                    d.crossVal = [randperm(size(d.label,1)/2) randperm(size(d.label,1)/2)]';
                case 'autoRun'
                    %Leave-one-run-pair-out, cond label balanced only: Random
                    %selection of N examples for the testing set at each
                    %fold, N being the number of observation per run times
                    %two grating conditions. Each example goes into the
                    %testing set once and only once. Does not care about
                    %sessions.
                    nSplit = str2double(p.split(1));
                    tmp1 = [];
                    tmp2 = [];
                    for i = 1:nSplit
                        tmp1 = [tmp1 randperm(size(d.label,1)/2/nSplit)];
                        tmp2 = [tmp2 randperm(size(d.label,1)/2/nSplit)];
                    end
                    tmp1 = tmp1(randperm(length(tmp1)));
                    tmp2 = tmp2(randperm(length(tmp2)));
                    d.crossVal = cat(1,tmp1',tmp2');
                case 'autoCmplt'
                    %Leave-one-pair-out, cond label balanced only:
                    %Deterministic selection of two consecutively acqired
                    %examples for the testing set at each fold. Each
                    %example goes into the testing set once and only once.
                    %Does not care about sessions.
                    d.crossVal = [1:size(d.label,1)/2 1:size(d.label,1)/2]';
                case 'sessBal'
                    %Leave-one-pair-out-per-session, cond label balanced:
                    %Random selection of four examples
                    %(sess1cond1+sess2cond1+sess1cond2+sess2cond2 + ... if more
                    %sessions) for the testing set at each fold. Each example
                    %goes into the testing set once and only once.
                    if strcmp(p.preSelection,'sess1') || strcmp(p.preSelection,'sess2')
                        error([p.preSelection ' does not make sense with ' p.k])
                    end
                    keyboard
                otherwise
                    error('')
            end
        else
            error('')
        end
    end
    
%     %the old way, counter-balancing a lot of shit
%     [dx,px] = handleKfolding(d,p);

%     %leave one pair out
%     d.crossVal = repmat(1:length(find(d.label==1)),[1 2])';

%     %leave one random pair out
%     d.crossVal = [randperm(length(find(d.label==1))) randperm(length(find(d.label==2)))]';
%         keyboard
      %purely random within runSVM3.m
%       d.crossVal = ([1:size(d.label)/2 1:size(d.label)/2])';
    

    %Run SVM
    
    switch p.algorithm
        case 'runSVM'
            switch p.keepInfo
                case 'p_then_m'
                    %Run on phase
                    pCur = p;
                    pCur.keepInfo = 'p';
                    [svm.r.hitRate(rep,:), ~,sortInd,svm.w(:,:,:,rep),svm.a(:,:,:,rep),svm.mi(:,:,rep),svm.data(:,:,:,rep)] = runSVM(d,pCur,verbose);
                    
                    %Run on magnitude with the same k-folding and sorting
                    %as on phase.
                    pCur = p;
                    pCur.keepInfo = 'm';
                    [svm.alt.r.hitRate(rep,:), ~,~,svm.alt.w(:,:,:,rep),svm.alt.a(:,:,:,rep),~,svm.alt.data(:,:,:,rep)] = runSVM(d,pCur,verbose,sortInd);
                case 'm_then_p'
                    %Run on mag
                    pCur = p;
                    pCur.keepInfo = 'm';
                    [svm.r.hitRate(rep,:), ~,sortInd,svm.w(:,:,:,rep),svm.a(:,:,:,rep),svm.mi(:,:,rep),svm.data(:,:,:,rep)] = runSVM(d,pCur,verbose);
                    
                    %Run on phase with the same k-folding and sorting
                    %as on mag.
                    pCur = p;
                    pCur.keepInfo = 'p';
                    [svm.alt.r.hitRate(rep,:), ~,~,svm.alt.w(:,:,:,rep),svm.alt.a(:,:,:,rep),~,svm.alt.data(:,:,:,rep)] = runSVM(d,pCur,verbose,sortInd);
                case {'p','m'}
%                     keyboard
%                     [svm.r.hitRate(rep,:), ~,sortInd,svm.w(:,:,:,rep),svm.a(:,:,:,rep),svm.mi(:,:,rep)] = runSVM(d,p,verbose);
%                     allRot = 0:10:180;
%                     for i = 1:length(allRot)
%                         p.rotTo = allRot(i);
%                         [curHitRate, ~,~,curw,cura,~] = runSVM(d,p,verbose,sortInd);
%                         allHitRate(i,:) = curHitRate;
%                         allW(i,:) = squeeze(nanmean(abs(curw),3));
%                         allA(i,:) = squeeze(nanmean(abs(cura),3));
%                     end
%                     allW
%                     figure('WindowStyle','docked')
%                     imagesc(allW)
%                     imagesc(allA)
                    
                    
                    if p.thresholdData
                        
                        [svm.r.hitRate(rep,:,:)                   , ~,~,~               ,~               ,~              ,~,tmpxValCorr          ] = runSVM4(d,p,verbose);
%                         [svm.r.hitRate                   , ~,~,~               ,~               ,~              ,~,tmpxValCorr          ] = runSVM4(d,p,verbose);
%                         tmp = squeeze(svm.r.hitRate);
%                         imagesc(tmp,[0 1]);
%                         set(gca,'YDir','reverse')
%                         set(gca,'XDir','reverse')
%                         xlabel('fun')
%                         ylabel('feat')
                        
                        
%                         [svm.r.hitRate(rep,:), ~,~,svm.w(:,:,:,rep),svm.a(:,:,:,rep),svm.mi(:,:,rep),~,tmpxValCorr] = runSVM3(d,p,verbose);
                        xValCorr.Rsquared(rep,:) = tmpxValCorr.Rsquared;
                        xValCorr.RHO(rep,:) = tmpxValCorr.RHO;
                    else
                        [svm.r.hitRate(rep,:), ~,~,svm.w(:,:,:,rep),svm.a(:,:,:,rep),~,~,tmpxValCorr] = runSVM4(d,p,verbose,p.funcROI.vec.sortInd);
                        xValCorr.Rsquared(rep,:) = tmpxValCorr.Rsquared;
                        xValCorr.RHO(rep,:) = tmpxValCorr.RHO;
                    end
                otherwise
                    error
            end
        case 'runSVMcomplexNregular'
            error('doubleCheck')
            [svm.r.hitRate(rep,:), ~,~] = runSVMcomplexNregular(d,p,verbose);
%         case {'runSVMcomplex'}
%             [svm.r.hitRate(rep,:), ~,~] = runSVMcomplex(d,p,verbose);
%         case 'runSVM'
%             [svm.r.hitRate(rep,:), ~,~] = runSVM(d,p,verbose);
        otherwise
            error
    end
    if ~isempty(err)
        if verbose
            save(fullfile(dataDirOut,dataFileOut,'__err'),'err')
            stopThis = 1;
            err = [];
            break
        else
            error(err)
        end
    end
    
    if verbose
        dur = toc;
        display(['Repetition ' num2str(rep) '; Took ' num2str(round(dur)) 'sec'])
    end
    
    
    
end
svm.p = p;

