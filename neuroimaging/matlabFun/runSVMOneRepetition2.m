function [rTr,rTe] = runSVMOneRepetition2(d,p,verbose,rep)

if ~exist('verbose','var')
    verbose = 0;
end
% if ~exist('doPerm','var')
%     doPerm = 0;
% end

%% K-folding and permutations
[d.crossVal,d.label] = defineKnPerm(d,p);
% tmp = unique(d.crossVal);
% tmp2 = [];
% for i = 1:length(tmp)
%     tmp2(i,:) = [length(find(tmp(i)==d.crossVal(1:end/2))) length(find(tmp(i)==d.crossVal(end/2+1:end)))];
% end

% %% Define k-folding
% if ischar(p.k)
%     switch p.k
%         %         case {'auto','autoCmplt'}
%         %             k = length(d.crossVal)/2;
% %         case 'autoRun'
% %             d.crossVal = defineCrossVal(d,p);
% %             k = length(d.crossVal)/2/str2double(p.split(1));
%         case {'run','runRand','sess','8FoldsRunRand','leaveTwoRunOut'}
%             d.crossVal = defineCrossVal(d,p);
%             k = length(unique(d.crossVal));
%         otherwise
%             error('');
%     end
% else
%     error('')
%     k = p.k;
%     p.k = 'randK';
% end
% 
% %% Do permutations depending on kFolding
% if p.doPerm
%     if ischar(p.k)
%         switch p.k
%             case {'runRand' '8FoldsRunRand','run'}
%                 
%                 %Shuffle labels on a run by run basis
%                 %02
% %                 tmpLabel = [];
% %                 for run = 1:length(d.label)/8/2
% %                     tmpLabel = cat(1,tmpLabel,repmat(randperm(2),8,1));
% %                 end
% %                 d.label = cat(1,tmpLabel(:,1),tmpLabel(:,2));
% 
%                 %03
%                 nTrial = str2double(p.split(1));
%                 if p.splitHalf
%                     nTrial = nTrial/2;
%                 end
%                 tmpLabel = cat(1,ones(length(d.label)/nTrial/2,1),ones(length(d.label)/nTrial/2,1)*2);
%                 tmpLabel = tmpLabel(randperm(length(tmpLabel)));
%                 d.label = reshape(repmat(tmpLabel,1,nTrial)',[1 length(d.label)])';
%                 
% %                 for curCrossVal = 1:k
% %                     d.label(d.crossVal==curCrossVal)'
% %                 end
% 
% %                 %Shuffle labels within each kFold
% %                 for curCrossVal = 1:k
% %                     tmp = d.label(d.crossVal==curCrossVal);
% %                     d.label(d.crossVal==curCrossVal) = tmp(randperm(length(tmp)));
% %                 end
%             case 'run'
%             otherwise
%                 error('did not code for that')
%         end
%     else
%         error('don''t know about that')
%     end
% end


% %% Permute labels if needed (permRun)
% if p.doPerm
%     %Contigent on crossVal
%     for i = 1:k
%         curInd = find(d.crossVal==i);
%         d.label(curInd) = d.label(curInd(randperm(length(curInd))));
%     end
% end

% %% Permute labels if needed (permSwap)
% if p.doPerm
%     for i = 1:k
%         curInd = find(d.crossVal==i);
%         if randi(2)-1
%             d.label(curInd) = -(d.label(curInd)-1.5)+1.5;
%         end
%     end
% end

% %% Permute labels if needed (permData)
% if p.doPerm
%     d.xData = d.xData(randperm(size(d.xData,1)),:);
% end

% %% Permute labels if needed (permLabel)
% if p.doPerm
%     d.label = d.label(randperm(size(d.label,1)),:);
% end


%% Multilevel svm
if ~iscell(p.infoComb)
    p.infoComb = {p.infoComb};
end
    
for infoCombInd = 1:length(p.infoComb)
    curP = p; curP.infoComb = p.infoComb{infoCombInd};
    switch curP.infoComb
%         case {'all','all_cross1','all_cross2','ALrev','ALrev_cross1','ALrev_cross2'}
%             tmpStr = strsplit(curP.infoComb,'_'); if strcmp(tmpStr{end}(1:5),'cross'); curP.subCrossVal = str2num(tmpStr{end}(6)); else curP.subCrossVal = 0; end
%             [rTr{infoCombInd},rTe{infoCombInd}] = runMultiLevelSVM2(d,curP);

        case {'ALrev1','ALrev2','all1','combInfo_ALrev'}
            tmpStr = strsplit(curP.infoComb,'_'); if ~isempty(strfind(tmpStr{end},'cross')); curP.subCrossVal = str2num(tmpStr{end}(6)); else curP.subCrossVal = 'nan'; end
%             [rTr{infoCombInd},rTe{infoCombInd}] = runMultiLevelSVM2(d,curP);
%             [rTr{infoCombInd},rTe{infoCombInd}] = runMultiLevelSVM4(d,curP);
            [rTr{infoCombInd},rTe{infoCombInd}] = runMultiLevelSVM5(d,curP);
%             [rTr{infoCombInd},rTe{infoCombInd}] = runMultiLevelSVM6(d,curP);
        case 'none'
            error('X')
        otherwise
            error('X')
    end
end



%% Lighten this up
rTr = reduceR(rTr,1);
% rTr = reduceR(rTr,2);
% rTr = reduceR(rTr,3);

rTe = reduceR(rTe,1);
% rTe = reduceR(rTe,2);
% rTe = reduceR(rTe,3);

% ind = 1
% length(find(rTe{ind}.yhat>0==[ones(length(rTe{ind}.yhat)/2,1); zeros(length(rTe{ind}.yhat)/2,1)]))/length(rTe{ind}.yhat)
% rTe{ind}.acc
% ind = 2
% length(find(rTe{ind}.yhat>0==[ones(length(rTe{ind}.yhat)/2,1); zeros(length(rTe{ind}.yhat)/2,1)]))/length(rTe{ind}.yhat)
% rTe{ind}.acc

% function [hitRate,decision,svmBias] = multiComb(data,label,crossVal,p)
% 
% k = length(unique(crossVal));
% decision = nan(p.nObs,1);
% hitRate = nan(k,1);
% svmBias = nan(p.nObs,1);
% for fold = 1:k
%     %Define train and test data
%     trInd = crossVal~=fold;
%     teInd = ~trInd;
%     
%     tr = data(trInd,:);
%     te = data(teInd,:);
%     trLabel = label(trInd,:);
%     teLabel = label(teInd,:);
%     
%     %Norm
%     trShift = mean(tr,1);
%     tr = tr-repmat(trShift,size(tr,1),1);
%     te = te-repmat(trShift,size(te,1),1);
%     trScale = std(tr,[],1);
%     tr = tr./repmat(trScale,size(tr,1),1);
%     te = te./repmat(trScale,size(te,1),1);
%     
%     %Train
%     svmStruct = svmtrain(trLabel,tr,'-s 0 -t 0 -c 1 -q');
%     
%     %Test
%     [~, tmpHitRate, curDecision,~] = svmpredict4(teLabel,te,svmStruct,tr,{},p);
%     
%     %Compile
%     decision(teInd,1) = curDecision;
%     svmBias(teInd,1) = svmStruct.rho;
%     hitRate(fold,1) = tmpHitRate/100;
% end
% 
% function [hitRate,foldDecision] = multiCombWithCrossVal(dOrig,p,k,rep,verbose)
% if p.doCrossInfo
%     error('not meant for that')
% end
% 
% 
% foldDecision = nan(p.nObs,length(p.trainInfo)+1);
% hitRate = nan(k,length(p.trainInfo)+1);
% 
% for fold = 1:k
%     %Make this fit the old scheme
%     dSubFold = dOrig;
%     dSubFold.xData(dSubFold.crossVal==fold,:) = [];
%     dSubFold.normData(dSubFold.crossVal(1:end/2)==fold,:) = []; % perpetuating that non-sense for backward compatibility
%     dSubFold.label(dSubFold.crossVal==fold) = [];
%     dSubFold.sessionLabel(dSubFold.crossVal==fold) = [];
%     dSubFold.crossVal(dSubFold.crossVal==fold) = [];
%     dSubFold.crossVal(dSubFold.crossVal>fold) = dSubFold.crossVal(dSubFold.crossVal>fold)-1;
%     pSubFold = p;
%     pSubFold.nObs = length(find(dOrig.crossVal~=fold));
%     
%     %Run the subLoop
%     %Get decision value fro train data
%     subFoldDecision = nan(pSubFold.nObs,length(p.trainInfo));
%     for subFold = 1:k-1
%         [~,curDecision,teInd,~,~,~] = trainNtest5(dSubFold,pSubFold,subFold,k-1,rep,verbose,p.trainInfo{1},p.testInfo{1});
%         subFoldDecision(teInd,1) = curDecision;
%         [~,curDecision,teInd,~,~,~] = trainNtest5(dSubFold,pSubFold,subFold,k-1,rep,verbose,p.trainInfo{2},p.testInfo{2});
%         subFoldDecision(teInd,2) = curDecision;
%     end
%     %Train SVM on all train data in from current fold (firtLevelSVMfilter)
%     %Apply firtLevelSVMfilter on test data from current fold (foldDecision)
%     [hitRate(fold,1),curDecision,teInd,~,~,~] = trainNtest5(dOrig,p,fold,k,rep,verbose,p.trainInfo{1},p.testInfo{1});
%     foldDecision(teInd,1) = curDecision;
%     [hitRate(fold,2),curDecision,teInd,~,~,~] = trainNtest5(dOrig,p,fold,k,rep,verbose,p.trainInfo{2},p.testInfo{2});
%     foldDecision(teInd,2) = curDecision;
%     %Train SVM on subFoldDecision in current fold (secondLevelSVMfilter)
%     %Apply secondLevelSVMfilter on foldDecision from current fold
%     tr = subFoldDecision;
%     te = foldDecision(teInd,1:length(p.trainInfo));
%     
%     %Norm
%     trShift = mean(tr,1);
%     tr = tr-repmat(trShift,size(tr,1),1);
%     te = te-repmat(trShift,size(te,1),1);
%     trScale = std(tr,[],1);
%     tr = tr./repmat(trScale,size(tr,1),1);
%     te = te./repmat(trScale,size(te,1),1);
%     
%     %Train
%     trLabel = ones(size(tr,1),1); trLabel(end/2+1:end) = 2;
%     teLabel = ones(size(te,1),1); teLabel(end/2+1:end) = 2;
%     svmStruct = svmtrain(trLabel,tr,'-s 0 -t 0 -c 1 -q');
%     
%     %Test
% %     [~, tmp, ~] = svmpredict(teLabel,te,svmStruct,'-q');
% %     hitRate(fold,1) = tmp(1);    
%     [~, tmpHitRate, curDecision,out] = svmpredict4(teLabel,te,svmStruct,tr,{},p);
%     foldDecision(teInd,length(p.trainInfo)+1) = curDecision;
%     hitRate(fold,length(p.trainInfo)+1) = tmpHitRate/100;
% end
% hitRate = mean(hitRate,1);
% 
% 
% function [hitRate,decision_values,teInd,filterdData,w,a,svmBias] = trainNtest5(dOrig,p,fold,k,rep,verbose,trainInfo,testInfo)
% %% Train
% % Prepare train data
% if ~strcmp(p.k,'sess')
%     [tr,~,teInd1,teInd2] = prepareDataAtFold(dOrig,p,trainInfo,fold,k,rep,verbose);
% else
%     [tr,te,teInd1,teInd2] = prepareDataAtFold(dOrig,p,trainInfo,fold,k,rep,verbose);
% end
% 
% % %Permute tr labels if needed (permFold)
% % if p.doPerm
% %     tr.label = tr.label(randperm(length(tr.label)));
% % end
% 
% 
% switch trainInfo
%     case {'m','mp','p','i','r'}
%         [X,Y] = pol2cart(tr.delay,tr.amp);
%         trFinal.data = [X Y]; clear X Y
%         trFinal.label = tr.label;
%     case 'm_i'
%         trFinal.data = [tr{1}.amp tr{2}.amp];
%         trFinal.label = tr{1}.label;
%     case 'delay'
%         trFinal.data = tr.delay;
%         trFinal.label = tr.label;
%     case {'amp','ampAtDelay','ampAtVoxDelay'}
%         trFinal.data = tr.amp;
%         trFinal.label = tr.label;
%     otherwise
%         error('XX')
% end
% 
% 
% 
% teInd = [teInd1; teInd2];
% 
% 
% % Train SVM
% if ~any(~strcmp(p.kernel,{'linear','rbf'}))
%     error('did not code for that')
% end
% % svmStruct = svmtrain(trFinal.data,trFinal.label,'kernel_function',p.kernel,'autoscale',0,'method','SMO','boxconstraint',p.C,'options',statset('MaxIter',1500000,'Display','off'));
% 
% 
% switch p.kernel
%     case 'linear'
%         kern = '0';
%     case 'rbf'
%         kern = '2';
%     otherwise
%         error('XX')
% end
% svmStruct = svmtrain(trFinal.label,trFinal.data,['-s 0 -t ' kern ' -c ' num2str(p.C) ' -q']);
% svmBias = svmStruct.rho;
% %% Test
% % Prepare test data
% if ~strcmp(p.k,'sess')
%     [~,te,~,~,~,tePartNorm] = prepareDataAtFold(dOrig,p,testInfo,fold,k,rep,verbose);
% end
% 
% switch trainInfo
%     case {'m','mp','p','i','r'}
%         [X,Y] = pol2cart(te.delay,te.amp);
%         teFinal.data = [X Y]; clear X Y
%         teFinal.label = te.label;
%     case 'm_i'
%         teFinal.data = [te{1}.amp te{2}.amp];
%         teFinal.label = te{1}.label;
%     case 'delay'
%         teFinal.data = te.delay;
%         teFinal.label = te.label;
%     case {'amp','ampAtDelay','ampAtVoxDelay'}
%         teFinal.data = te.amp;
%         teFinal.label = te.label;
%     otherwise
%         error('XX')
% end
% 
% 
% % if strcmp(p.kernel,'linear')
% %     [~, hitRate, decision_values,out] = svmpredict4(teFinal.label,teFinal.data,svmStruct,trFinal.data,{trainInfo testInfo},p);
% %     hitRate = hitRate/100;
% %     
% %     a = out.a;
% %     w = out.w;
% % else
%     [~, hitRate, decision_values] = svmpredict(teFinal.label,teFinal.data,svmStruct,'-q');
%     hitRate = hitRate(1)/100;
%     a = [];
%     w = [];
% % end
% 
% 
% % tePartNorm = finalPrep(tePartNorm,'mp');
% % filterdData = wFilter(tePartNorm.data,svmStruct,trFinal.data);
% filterdData = [];
% 
% 
% 
% 
% 
% function [w,a] = trainOnly(dOrig,p,fold,k,rep,verbose,trainInfo,testInfo)
% %% Train
% % Prepare train data
% dOrig.crossVal = zeros(size(dOrig.crossVal));
% if ~strcmp(p.k,'sess')
%     [tr,~,teInd1,teInd2] = prepareDataAtFold(dOrig,p,trainInfo,fold,k,rep,verbose);
% else
%     [tr,te,teInd1,teInd2] = prepareDataAtFold(dOrig,p,trainInfo,fold,k,rep,verbose);
% end
% 
% % %Permute tr labels if needed (permFold)
% % if p.doPerm
% %     tr.label = tr.label(randperm(length(tr.label)));
% % end
% 
% 
% switch trainInfo
%     case {'m','mp','p','i','r'}
%         [X,Y] = pol2cart(tr.delay,tr.amp);
%         trFinal.data = [X Y]; clear X Y
%         trFinal.label = tr.label;
%     case 'm_i'
%         trFinal.data = [tr{1}.amp tr{2}.amp];
%         trFinal.label = tr{1}.label;
%     case 'delay'
%         trFinal.data = tr.delay;
%         trFinal.label = tr.label;
%     case {'amp','ampAtDelay','ampAtVoxDelay'}
%         trFinal.data = tr.amp;
%         trFinal.label = tr.label;
%     otherwise
%         error('XX')
% end
% 
% 
% 
% % teInd = [teInd1; teInd2];
% 
% 
% % Train SVM
% if ~any(~strcmp(p.kernel,{'linear','rbf'}))
%     error('did not code for that')
% end
% % svmStruct = svmtrain(trFinal.data,trFinal.label,'kernel_function',p.kernel,'autoscale',0,'method','SMO','boxconstraint',p.C,'options',statset('MaxIter',1500000,'Display','off'));
% 
% 
% switch p.kernel
%     case 'linear'
%         kern = '0';
%     case 'rbf'
%         kern = '2';
%     otherwise
%         error('XX')
% end
% svmStruct = svmtrain(trFinal.label,trFinal.data,['-s 0 -t ' kern ' -c ' num2str(p.C) ' -q']);
% svmBias = svmStruct.rho;
% 
% [w,a] = getAandW(svmStruct,trFinal.data);
% 
% 
% 
% 
% function [tr,te,teInd1,teInd2,trPartNorm,tePartNorm] = prepareDataAtFold(d,p,keepInfo,fold,k,rep,verbose)
% 
% 
% if verbose
%     tic
%     display(['Fold ' num2str(fold) '/' num2str(k) '; Rep ' num2str(rep) '/' num2str(p.repeat)])
% end
% 
% %Initiate some stuff
% 
% %     wFold = w(:,:,:,fold);
% %     Afold = A(:,:,:,fold);
% 
% %Break complex data down to amp and delay
% d2.delay(:,:,1) = angle(d.xData(1:end/2,:));
% d2.delay(:,:,2) = angle(d.xData(end/2+1:end,:));
% d2.delay(:,:,3) = angle(d.normData);
% d2.amp(:,:,1) = abs(d.xData(1:end/2,:));
% d2.amp(:,:,2) = abs(d.xData(end/2+1:end,:));
% d2.amp(:,:,3) = abs(d.normData);
% % switch keepInfo
% %     case {'m','p','mp'}
% %     case {'r','i'}
% %         %remove voxel mean phase estimated from the plaid
% %         d2.delay(:,:,1) = wrapToPi(d2.delay(:,:,1)-repmat(circ_mean(d2.delay(:,:,3),[],1),[size(d2.delay,1) 1]));
% %         d2.delay(:,:,2) = wrapToPi(d2.delay(:,:,2)-repmat(circ_mean(d2.delay(:,:,3),[],1),[size(d2.delay,1) 1]));
% %         d2.delay(:,:,3) = wrapToPi(d2.delay(:,:,3)-repmat(circ_mean(d2.delay(:,:,3),[],1),[size(d2.delay,1) 1]));
% %         switch keepInfo
% %             case 'r'
% %                 %set amplitude to the real value (amplitude at the mean delay)
% %                 [d2.amp,~] = pol2cart(d2.delay,d2.amp);
% %             case 'i'
% %                 %set amplitude to the mag-normalized imaginary value (delay relative to mean delay)
% %                 [~,d2.amp] = pol2cart(d2.delay,1);
% %         end
% %     otherwise
% %         error('double-check that')
% % end
% d2.sessionLabel(:,:,1) = d.sessionLabel(1:end/2,:);
% d2.sessionLabel(:,:,2) = d.sessionLabel(end/2+1:end,:);
% d2.sessionLabel(:,:,3) = d.sessionLabel(1:end/2,:);
% d2.label(:,:,1) = d.label(1:end/2,:);
% d2.label(:,:,2) = d.label(end/2+1:end,:);
% d2.label(:,:,3) = ones(size(d.label(1:end/2,:))).*3;
% d2.crossVal(:,:,1) = d.crossVal(1:end/2,:);
% d2.crossVal(:,:,2) = d.crossVal(end/2+1:end,:);
% d2.crossVal(:,:,3) = nan(size(d.crossVal(1:end/2,:)));
% 
% %Split train and test set
% teInd1 = d2.crossVal(:,1,1)==fold;
% teInd2 = d2.crossVal(:,1,2)==fold;
% allFields = fields(d2);
% for i = 1:length(allFields)
%     d2te.(allFields{i})(:,:,1) = d2.(allFields{i})(teInd1,:,1);
%     d2te.(allFields{i})(:,:,2) = d2.(allFields{i})(teInd2,:,2);
%     tmp1 = d2.(allFields{i})(~teInd1,:,1);
%     tmp2 = d2.(allFields{i})(~teInd2,:,2);
%     d2tr.(allFields{i}) = cat(3,tmp1,tmp2); clear tmp1 tmp2
% end
% 
% 
% 
% % if p.regSession
% %     error('double-check that')
% %     %estimate session effect from training set
% %     [betaHat,X] = estimateSessionEffect(d2tr.amp,d2tr.sessionLabel);
% %     firsSessInd = size(X,2)-(length(unique(d2tr.sessionLabel))-1)+1;
% %     sessionEffect = betaHat(:,firsSessInd:end)';
% %     
% %     %remove session effect
% %     for sess = 2:length(unique(d2tr.sessionLabel))
% %         %from training set
% %         sessInd = d2tr.sessionLabel(:,:,1)==sess;
% %         d2tr.amp(sessInd,:,:) = d2tr.amp(sessInd,:,:)-repmat(sessionEffect,[length(find(sessInd)) 1 3]);
% %         %from testing set
% %         sessInd = d2te.sessionLabel(:,:,1)==sess;
% %         d2te.amp(sessInd,:,:) = d2te.amp(sessInd,:,:)-repmat(sessionEffect,[length(find(sessInd)) 1 3]);
% %     end
% % end
% 
% 
% 
% %% Normalize
% allFields = fields(d2tr);
% for i = 1:length(allFields)
%     tr.(allFields{i}) = cat(1,d2tr.(allFields{i})(:,:,1),d2tr.(allFields{i})(:,:,2));
% end
% allFields = fields(d2te);
% for i = 1:length(allFields)
%     te.(allFields{i}) = cat(1,d2te.(allFields{i})(:,:,1),d2te.(allFields{i})(:,:,2));
% end
% 
% % figure('WindowStyle','docked');
% % [X,Y] = pol2cart(tr.delay,tr.amp);
% % scatter(mean(X,1),mean(Y,1))
% 
% %Subtract mean delay voxel-wise
% if ~strcmp(keepInfo,'ampAtDelay') && ~strcmp(keepInfo,'ampAtVoxDelay')
%     phase_shift1 = circ_mean(tr.delay(1:end/2,:),[],1);
%     phase_shift2 = circ_mean(tr.delay(end/2+1:end,:),[],1);
%     phase_shift = circ_mean([phase_shift1; phase_shift2],[],1);
%     tr.delay = wrapToPi(tr.delay-repmat(phase_shift,size(tr.delay,1),1));
%     if strcmp(p.k,'sess')
%         error('double-check that')
%         phase_shift1 = circ_mean(te.delay(1:end/2,:),[],1);
%         phase_shift2 = circ_mean(te.delay(end/2+1:end,:),[],1);
%         phase_shift = circ_mean([phase_shift1; phase_shift2],[],1);
%     end
%     te.delay = wrapToPi(te.delay-repmat(phase_shift,size(te.delay,1),1));
%     %Scale amp to mean of one
%     amp_scale = mean(tr.amp,1);
%     tr.amp = tr.amp./repmat(amp_scale,[size(tr.amp,1) 1]);
%     if strcmp(p.k,'sess')
%         amp_scale = mean(te.amp,1);
%     end
%     te.amp = te.amp./repmat(amp_scale,[size(te.amp,1) 1]);
% end
% 
% % figure('WindowStyle','docked');
% % [X,Y] = pol2cart(tr.delay,tr.amp);
% % scatter(mean(X,1),mean(Y,1))
% 
% trPartNorm = tr;
% tePartNorm = te;
% 
% switch keepInfo
%     case 'r'
%         %Take real part
%         [tr.amp,~] = pol2cart(tr.delay,tr.amp);
%         [te.amp,~] = pol2cart(te.delay,te.amp);
%         %Set delay to 0
%         tr.delay = zeros(size(tr.delay));
%         te.delay = zeros(size(te.delay));
% %         %Z-score
% %         amp_shift = mean(tr.amp,1);
% %         tr.amp = tr.amp-repmat(amp_shift,size(tr.amp,1),1);
% %         te.amp = te.amp-repmat(amp_shift,size(te.amp,1),1);
% %         amp_scale = std(tr.amp,[],1);
% %         tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
% %         te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
% %     case 'i'
% %         %Take imaginary part
% %         [~,tr.amp] = pol2cart(tr.delay,tr.amp);
% %         [~,te.amp] = pol2cart(te.delay,te.amp);
% %         %Set delay to pi/2
% %         tr.delay = ones(size(tr.delay))*pi/2;
% %         te.delay = ones(size(te.delay))*pi/2;
% %         %         %Z-score
% %         %         amp_shift = mean(tr.amp,1);
% %         %         tr.amp = tr.amp-repmat(amp_shift,size(tr.amp,1),1);
% %         %         te.amp = te.amp-repmat(amp_shift,size(te.amp,1),1);
% %         %         amp_scale = std(tr.amp,[],1);
% %         %         tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
% %         %         te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
%     case 'i'
%         %Remove amp
%         tr.amp = ones(size(tr.amp));
%         te.amp = ones(size(te.amp));
%         %Take imaginary part
%         [~,tr.amp] = pol2cart(tr.delay,tr.amp);
%         [~,te.amp] = pol2cart(te.delay,te.amp);
%         %Set delay to pi/2
%         tr.delay = ones(size(tr.delay))*pi/2;
%         te.delay = ones(size(te.delay))*pi/2;
%         %         %Z-score
%         %         amp_shift = mean(tr.amp,1);
%         %         tr.amp = tr.amp-repmat(amp_shift,size(tr.amp,1),1);
%         %         te.amp = te.amp-repmat(amp_shift,size(te.amp,1),1);
%         %         amp_scale = std(tr.amp,[],1);
%         %         tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
%         %         te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
%     case {'m','amp'}
%         %Remove delay
%         tr.delay = zeros(size(tr.delay));
%         te.delay = zeros(size(te.delay));
% %         %Z-score)
% %         amp_shift = mean(tr.amp,1);
% %         tr.amp = tr.amp-repmat(amp_shift,size(tr.amp,1),1);
% %         te.amp = te.amp-repmat(amp_shift,size(te.amp,1),1);
% %         amp_scale = std(tr.amp,[],1);
% %         tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
% %         te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
%     case 'm_i'
%         trOrig = tr;
%         teOrig = te;
% 
%         %Do m
%         tr = trOrig;
%         te = teOrig;
%         %Remove delay
%         tr.delay = zeros(size(tr.delay));
%         te.delay = zeros(size(te.delay));
%         
%         tr_m = tr;
%         te_m = te;
%         
%         
%         %Do i
%         tr = trOrig;
%         te = teOrig;
%         %Remove amp
%         tr.amp = ones(size(tr.amp));
%         te.amp = ones(size(te.amp));
%         %Take imaginary part
%         [~,tr.amp] = pol2cart(tr.delay,tr.amp);
%         [~,te.amp] = pol2cart(te.delay,te.amp);
%         %Set delay to pi/2
%         tr.delay = ones(size(tr.delay))*pi/2;
%         te.delay = ones(size(te.delay))*pi/2;
%         
%         tr_i = tr;
%         te_i = te;
%         
%         
%         tr = trOrig;
%         te = teOrig;
%     case {'p','delay'}
%         %Remove amp
%         tr.amp = ones(size(tr.amp));
%         te.amp = ones(size(te.amp));
%     case 'mp'
% %         %Amp scale only (equivalent to image scaling only)
% %         amp_scale = mean(tr.amp,1);
% %         tr.amp = tr.amp./repmat(amp_scale,[size(tr.amp,1) 1]);
% %         te.amp = te.amp./repmat(amp_scale,[size(te.amp,1) 1]);
% %         %Image Shift and Scale
% %         [trR,trI] = pol2cart(tr.delay,tr.amp);
% %         [teR,teI] = pol2cart(te.delay,te.amp);
% %         %Shift
% %         %real
% %         trR_shift = mean(trR,1);
% %         trR = trR-repmat(trR_shift,size(trR,1),1);
% %         teR = teR-repmat(trR_shift,size(teR,1),1);
% %         %imag
% %         trI_shift = mean(trI,1);
% %         trI = trI-repmat(trI_shift,size(trI,1),1);
% %         teI = teI-repmat(trI_shift,size(teI,1),1);
% %         %Scale
% %         tr_scale = mean([std(trR,[],1); std(trI,[],1)],1);
% %         trR = trR./repmat(tr_scale,size(trR,1),1); trI = trI./repmat(tr_scale,size(trI,1),1);
% %         teR = teR./repmat(tr_scale,size(teR,1),1); teI = teI./repmat(tr_scale,size(teI,1),1);
% %         
% %         [tr.delay,tr.amp] = cart2pol(trR,trI);
% %         [te.delay,te.amp] = cart2pol(teR,teI);
% %         clear trR trI tr_shift tr_scale teR teI te_shift te_scale
%     case {'ampAtDelay','ampAtVoxDelay'}
%         
%     otherwise
%         error('did not code that')
% end
% 
% switch keepInfo
%     case {'m','mp','i','p','r'}
%         [tr,te] = doZ(tr,te,p);
%     case 'm_i'
%         [tr_i,te_i] = doZ(tr_i,te_i,p);
%         [tr_m,te_m] = doZ(tr_m,te_m,p);
%         tr = {tr_i tr_m};
%         te = {te_i te_m};
%     case 'delay'
%         %shit already performed earlier
%         delay_scale = circ_std(tr.delay,[],[],1);
%         tr.delay = tr.delay./repmat(delay_scale,size(tr.delay,1),1);
%         te.delay = te.delay./repmat(delay_scale,size(te.delay,1),1);
%     case 'amp'
%         amp_shift = mean(tr.amp,1);
%         tr.amp = tr.amp - repmat(amp_shift,size(tr.amp,1),1);
%         te.amp = te.amp - repmat(amp_shift,size(te.amp,1),1);
%         amp_scale = std(tr.amp,[],1);
%         tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
%         te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
%     case 'ampAtDelay'
% %         figure('WindowStyle','docked');
% %         [X,Y] = pol2cart(tr.delay,tr.amp);
% %         scatter(mean(X,1),mean(Y,1))
%         
% %         delay_shift1 = circ_mean(tr.delay(1:end/2,:),[],1);
% %         delay_shift2 = circ_mean(tr.delay(end/2+1:end,:),[],1);
% %         delay_shift = circ_mean([delay_shift1; delay_shift2],[],1);
% %         delay_shift = circ_median(delay_shift,2);
% %         tr.delay2 = wrapToPi(tr.delay - repmat(delay_shift,size(tr.delay,1),size(tr.delay,2)));
% %         figure('WindowStyle','docked');
% %         [X,Y] = pol2cart(tr.delay2,tr.amp);
% %         scatter(mean(X,1),mean(Y,1))
%         
%         % Project data at ROI median delay
%         delay_shift1 = circ_mean(tr.delay(1:end/2,:),[],1);
%         delay_shift2 = circ_mean(tr.delay(end/2+1:end,:),[],1);
%         delay_shift = circ_mean([delay_shift1; delay_shift2],[],1);
%         delay_shift = circ_median(delay_shift,2);
%         [tr.amp,~] = pol2cart(wrapToPi(tr.delay - repmat(delay_shift,size(tr.delay,1),size(tr.delay,2))),tr.amp);
% %         [~,tr.amp] = pol2cart(wrapToPi(tr.delay - repmat(delay_shift,size(tr.delay,1),size(tr.delay,2))),tr.amp);
%         tr.delay = zeros(size(tr.amp));
%         [te.amp,~] = pol2cart(wrapToPi(te.delay - repmat(delay_shift,size(te.delay,1),size(te.delay,2))),te.amp);
% %         [~,te.amp] = pol2cart(wrapToPi(te.delay - repmat(delay_shift,size(te.delay,1),size(te.delay,2))),te.amp);
%         te.delay = zeros(size(te.amp));
%         
% %         figure('WindowStyle','docked');
% %         [X,Y] = pol2cart(tr.delay,tr.amp);
% %         scatter(mean(X,1),mean(Y,1))
%         
%         % z-score
%         amp_shift = mean(tr.amp,1);
%         tr.amp = tr.amp - repmat(amp_shift,size(tr.amp,1),1);
%         te.amp = te.amp - repmat(amp_shift,size(te.amp,1),1);
%         amp_scale = std(tr.amp,[],1);
%         tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
%         te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
%         
% %         figure('WindowStyle','docked');
% %         [X,Y] = pol2cart(tr.delay,tr.amp);
% %         scatter(mean(X,1),mean(Y,1))
%         
%     case 'ampAtVoxDelay'
% %         figure('WindowStyle','docked');
% %         [X,Y] = pol2cart(tr.delay,tr.amp);
% %         scatter(mean(X,1),mean(Y,1))
% %         
% %         delay_shift1 = circ_mean(tr.delay(1:end/2,:),[],1);
% %         delay_shift2 = circ_mean(tr.delay(end/2+1:end,:),[],1);
% %         delay_shift = circ_mean([delay_shift1; delay_shift2],[],1);
% %         tr.delay = wrapToPi(tr.delay - repmat(delay_shift,size(tr.delay,1),1));
% %         figure('WindowStyle','docked');
% %         [X,Y] = pol2cart(tr.delay,tr.amp);
% %         scatter(mean(X,1),mean(Y,1))
%         
%         % Project data at voxel-specific delay
%         delay_shift1 = circ_mean(tr.delay(1:end/2,:),[],1);
%         delay_shift2 = circ_mean(tr.delay(end/2+1:end,:),[],1);
%         delay_shift = circ_mean([delay_shift1; delay_shift2],[],1);
%         [tr.amp,~] = pol2cart(wrapToPi(tr.delay - repmat(delay_shift,size(tr.delay,1),1)),tr.amp);
% %         [~,tr.amp] = pol2cart(wrapToPi(tr.delay - repmat(delay_shift,size(tr.delay,1),1)),tr.amp);
%         tr.delay = zeros(size(tr.amp));
%         [te.amp,~] = pol2cart(wrapToPi(te.delay - repmat(delay_shift,size(te.delay,1),1)),te.amp);
% %         [~,te.amp] = pol2cart(wrapToPi(te.delay - repmat(delay_shift,size(te.delay,1),1)),te.amp);
%         te.delay = zeros(size(te.amp));
%         % z-score
%         amp_shift = mean(tr.amp,1);
%         tr.amp = tr.amp - repmat(amp_shift,size(tr.amp,1),1);
%         te.amp = te.amp - repmat(amp_shift,size(te.amp,1),1);
%         amp_scale = std(tr.amp,[],1);
%         tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
%         te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
%     otherwise
%         error('XX')
% end
% 
% 
% 
% 
% 
% function xData = finalPrep(data,keepInfo)
% % switch keepInfo
% %     case {'m','r','i','iN'}
% %         xData.data = data.amp;
% %         xData.label = data.label;
% %     case {'p','mp'}
% %         [X,Y] = pol2cart(data.delay,data.amp);
% %         xData.data = [X Y];
% %         xData.label = data.label;
% % %     case 'mp'
% % %         [X,Y] = pol2cart(data.delay,data.amp);
% % %         xData.data = [X Y];
% % %         xData.label = data.label;
% %     otherwise
% %         error('xx')
% % end
% [X,Y] = pol2cart(data.delay,data.amp);
% xData.data = [X Y];
% xData.label = data.label;
% 
% 
% 
% function [predictedLabel,accuracy,d,out] = svmpredict4(y,x,svmStruct,xTr,info,p)
% 
% [w,a] = getAandW(svmStruct,xTr);
% 
% 
% %% Adapt w to test data if needed
% % keyboard
% % 
% % wC = complex(w(1:end/2),w(end/2+1:end));
% % aC = complex(a(1:end/2),a(end/2+1:end));
% % 
% % close all
% % figure('WindowStyle','docked');
% % h = scatter(angle(wC),abs(wC),'filled','k'); hold on
% % alpha(h,.2);
% % plot([-pi/2 -pi/2],get(gca,'ylim'),'r')
% % plot([pi/2 pi/2],get(gca,'ylim'),'r')
% % title('w')
% % 
% % figure('WindowStyle','docked');
% % h = scatter(angle(aC),abs(aC),'filled','k'); hold on
% % alpha(h,.2);
% % plot([-pi/2 -pi/2],get(gca,'ylim'),'r')
% % plot([pi/2 pi/2],get(gca,'ylim'),'r')
% % title('a')
% % 
% % figure('WindowStyle','docked');
% % hist([angle(wC); angle(aC)]',100)
% % legend({'w','a'})
% % 
% % 
% % lims = [min([real(wC) imag(wC)]) max([real(wC) imag(wC)])];
% % h = scatter(real(wC),imag(wC),'filled','k'); hold on
% % alpha(h,.2);
% % xlim(lims); ylim(lims);
% % plot(get(gca,'xlim'),[0 0],'k')
% % plot([0 0],get(gca,'ylim'),'k')
% % 
% % 
% % 
% % figure('WindowStyle','docked');
% % lims = [min([real(wC) imag(wC)]) max([real(wC) imag(wC)])];
% % h = scatter(real(wC),imag(wC),'filled','k'); hold on
% % alpha(h,.2);
% % xlim(lims); ylim(lims);
% % plot(get(gca,'xlim'),[0 0],'k')
% % plot([0 0],get(gca,'ylim'),'k')
% 
% if p.doCrossInfo
%     w = crossInfoAdapt(w,info,p);
% end
% 
% % if size(x,2)==size(svmStruct.SVs,2)
% % elseif size(svmStruct.SVs,2)==size(x,2)*2
% %     w = sqrt(w(1,1:end/2).^2 + w(1,end/2+1:end).^2);
% %     w1 = sqrt(w1(1,1:end/2).^2 + w1(1,end/2+1:end).^2);
% %     w2 = sqrt(w2(1,1:end/2).^2 + w2(1,end/2+1:end).^2);
% % elseif size(svmStruct.SVs,2)==size(x,2)/2
% %     w = repmat(w,[1 2]);
% %     w1 = repmat(w1,[1 2]);
% %     w2 = repmat(w2,[1 2]);
% % else
% %     error('xx')
% % end
% 
% %% Get decision value d
% % d
% d = nan(size(x,1),1);
% for j=1:size(x,1)
%     d(j,1) = dot(w,x(j,:)) - svmStruct.rho;
% end
% % % d1
% % d1 = nan(size(x,1),1);
% % for j=1:size(x,1)
% %     d1(j,1) = dot(w1,x(j,:)) - svmStruct.rho/2;
% % end
% % % d2
% % d2 = nan(size(x,1),1);
% % for j=1:size(x,1)
% %     d2(j,1) = dot(w2,x(j,:)) - svmStruct.rho/2;
% % end
% 
% %% Predict label
% predictedLabel = ones(size(d));
% predictedLabel(d<0) = 2;
% 
% %% Compute accuracy
% accuracy = length(find(y==predictedLabel))/length(y)*100;
% 
% %% Compile out
% out.w = w;
% % out.w1 = w1;
% % out.w2 = w2;
% out.a = a;
% % out.a1 = a1;
% % out.a2 = a2;
% out.d = d;
% % out.d1 = d1;
% % out.d2 = d2;
% 
% 
% function [w,a] = getAandW(svmStruct,xTr)
% nTr = size(xTr,1);
% 
% %% Get w
% % w
% w = zeros(1,size(svmStruct.SVs,2));
% for i = 1:size(svmStruct.SVs,1)
%     w = w + svmStruct.sv_coef(i)*svmStruct.SVs(i,:);
% end
% % % w1
% % ind = find(svmStruct.sv_indices<=nTr/2);
% % w1 = zeros(1,size(svmStruct.SVs,2));
% % for i = 1:length(ind)% size(svmStruct.SVs,1)/2
% %     w1 = w1 + svmStruct.sv_coef(ind(i))*svmStruct.SVs(ind(i),:);
% % end
% % % w2
% % ind = find(svmStruct.sv_indices>nTr/2);
% % w2 = zeros(1,size(svmStruct.SVs,2));
% % for i = 1:length(ind)%size(svmStruct.SVs,1)/2+1:size(svmStruct.SVs,1)
% %     w2 = w2 + svmStruct.sv_coef(ind(i))*svmStruct.SVs(ind(i),:);
% % end
% 
% %% Get a
% % a
% s = w*xTr';
% % a = w*cov(xTr)/cov(s);
% a = nan(size(w));
% for dim = 1:size(xTr,2)
%     tmp = cov(xTr(:,dim),s);
%     a(dim) = tmp(1,2);
% end
% 
% % % a1
% % s1 = w1*xTr(1:end/2,:)';
% % % a1 = w1*cov(xTr(1:end/2,:))/cov(s1);
% % a1 = nan(size(w1));
% % for dim = 1:length(w2)
% %     tmp = cov(xTr(1:end/2,dim),s1);
% %     a1(dim) = tmp(1,2);
% % end
% % 
% % % a2
% % s2 = w2*xTr(end/2+1:end,:)';
% % % a2 = w2*cov(xTr(end/2+1:end,:))/cov(s2);
% % a2 = nan(size(w2));
% % for dim = 1:length(w2)
% %     tmp = cov(xTr(end/2+1:end,dim),s2);
% %     a2(dim) = tmp(1,2);
% % end
% 
% 
% function filteredData = wFilter(x,svmStruct,xTr)
% 
% nTr = size(xTr,1);
% 
% %% Get w
% % w
% w = zeros(1,size(svmStruct.SVs,2));
% for i = 1:size(svmStruct.SVs,1)
%     w = w + svmStruct.sv_coef(i)*svmStruct.SVs(i,:);
% end
% % w1
% ind = find(svmStruct.sv_indices<=nTr/2);
% w1 = zeros(1,size(svmStruct.SVs,2));
% for i = 1:length(ind)% size(svmStruct.SVs,1)/2
%     w1 = w1 + svmStruct.sv_coef(ind(i))*svmStruct.SVs(ind(i),:);
% end
% % w2
% ind = find(svmStruct.sv_indices>nTr/2);
% w2 = zeros(1,size(svmStruct.SVs,2));
% for i = 1:length(ind)%size(svmStruct.SVs,1)/2+1:size(svmStruct.SVs,1)
%     w2 = w2 + svmStruct.sv_coef(ind(i))*svmStruct.SVs(ind(i),:);
% end
% 
% 
% %% Split w and data in real and imaginary parts
% if size(w,2)==size(x,2)
%     wX = w(1:end/2);
%     wY = w(end/2+1:end);
%     w1X = w1(1:end/2);
%     w1Y = w1(end/2+1:end);
%     w2X = w2(1:end/2);
%     w2Y = w2(end/2+1:end);
% elseif size(w,2)==size(x,2)/2
%     wX = w;
%     wY = w;
%     w1X = w1;
%     w1Y = w1;
%     w2X = w2;
%     w2Y = w2;
% else
%     error('xx')
% end
% xX = x(:,1:end/2);
% xY = x(:,end/2+1:end);
% 
% %% Apply w filter
% % d
% fX = nan(size(xX,1),1);
% for j=1:size(xX,1)
%     fX(j,1) = dot(wX,xX(j,:));
% end
% fY = nan(size(xY,1),1);
% for j=1:size(xY,1)
%     fY(j,1) = dot(wY,xY(j,:));
% end
% filteredData.f = complex(fX,fY);
% % d1
% f1X = nan(size(xX,1),1);
% for j=1:size(xX,1)
%     f1X(j,1) = dot(w1X,xX(j,:));
% end
% f1Y = nan(size(xY,1),1);
% for j=1:size(xY,1)
%     f1Y(j,1) = dot(w1Y,xY(j,:));
% end
% filteredData.f1 = complex(f1X,f1Y);
% % d2
% f2X = nan(size(xX,1),1);
% for j=1:size(xX,1)
%     f2X(j,1) = dot(w2X,xX(j,:));
% end
% f2Y = nan(size(xY,1),1);
% for j=1:size(xY,1)
%     f2Y(j,1) = dot(w2Y,xY(j,:));
% end
% filteredData.f2 = complex(f2X,f2Y);

function crossVal = defineCrossVal(d,p)
        
switch p.k
    case 'run'
        if p.splitHalf
            nSplit = str2double(p.split(1))/2;
        else
            nSplit = str2double(p.split(1));
        end
        tmp = repmat(1:p.nObs/nSplit/2,nSplit,1);
        tmp1 = tmp(:,randperm(size(tmp,2)));
        tmp2 = tmp(:,randperm(size(tmp,2)));
        crossVal = cat(1,reshape(tmp1,p.nObs/2,1),reshape(tmp2,p.nObs/2,1));
    case 'leaveTwoRunOut'
        nSplit = str2double(p.split(1));
        tmp = 1:p.nObs/nSplit/2;
        tmp1 = tmp(randperm(length(tmp)));
        tmp2 = tmp(randperm(length(tmp)));
        tmp1 = repmat(tmp1,nSplit,1);
        tmp2 = repmat(tmp2,nSplit,1);
        tmp1 = reshape(tmp1,p.nObs/2,1);
        tmp2 = reshape(tmp2,p.nObs/2,1);
        crossVal = [tmp1; tmp2]; clear tmp1 tmp2
    case 'runRand'
        nSplit = str2double(p.split(1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nSplit = nSplit/2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp1 = [];
        tmp2 = [];
        for i = 1:nSplit
            tmp1 = [tmp1 randperm(size(d.label,1)/2/nSplit)];
            tmp2 = [tmp2 randperm(size(d.label,1)/2/nSplit)];
        end
        tmp1 = tmp1(randperm(length(tmp1)));
        tmp2 = tmp2(randperm(length(tmp2)));
        crossVal = cat(1,tmp1',tmp2');
    case '8FoldsRunRand'
        tmp1 = [];
        tmp2 = [];
        for i = 1:(size(d.label,1)/2/8)
            tmp1 = [tmp1 randperm(8)];
            tmp2 = [tmp2 randperm(8)];
        end
        crossVal = cat(1,tmp1',tmp2');
    case 'sess'
        crossVal = d.sessionLabel;        
end
