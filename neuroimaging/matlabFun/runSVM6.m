function hitRate = runSVM6(d,p,verbose,rep,doPerm)
if ~exist('verbose','var')
    verbose = 0;
end
if ~isfield(p,'rotTo')
    if strcmp(p.keepInfo,'i')
        p.rotTo = 0;
    else
        p.rotTo = pi/4;
    end
end
if ~exist('doPerm','var')
    doPerm = 0;
end

%% Permute data if needed
if doPerm
    d.xData = d.xData(randperm(size(d.xData,1)),:);
end

%% Define k-folding
if ischar(p.k)
    switch p.k
        %         case {'auto','autoCmplt'}
        %             k = length(d.crossVal)/2;
        case 'autoRun'
            d.crossVal = defineCrossVal(d,p);
            k = length(d.crossVal)/2/str2double(p.split(1));
        otherwise
            error('');
    end
else
    error('')
    k = p.k;
    p.k = 'randK';
end


%% Loop over cross-validation folds
dOrig = d;
switch p.algorithm
    case 'runSVM'
        hitRate = nan(k,length(p.featLevList),length(p.featLevList));
    case 'runSVM_RFE'
%         hitRate = nan(k,p.RFE_R,length(p.featLevList),p.RFE_L);
        hitRate = nan(k,p.RFE_R,length(p.featLevList));
%         w = nan(k,p.RFE_R,length(p.featLevList),p.RFE_L,p.nFeatures);
    otherwise
        error('did not code for that')
end
% w = nan(length(p.featLevList),length(p.featLevList),p.featLevList(end),k);
% A = nan(length(p.featLevList),length(p.featLevList),p.featLevList(end),k);
% % data = nan(length(p.featLevList),p.featLevList(end),k);
% sortInd = nan(k,p.nFeatures);
% mi = nan(k,p.nFeatures);
for fold = 1:k
    if verbose
        tic
        display(['Fold ' num2str(fold) '/' num2str(k) '; Rep ' num2str(rep) '/' num2str(p.repeat)])
    end
    % Initiate some stuff
    d = dOrig;
    
    %     wFold = w(:,:,:,fold);
    %     Afold = A(:,:,:,fold);
    
    % Break complex data down to amp and delay
    d2.delay(:,:,1) = angle(d.xData(1:end/2,:));
    d2.delay(:,:,2) = angle(d.xData(end/2+1:end,:));
    d2.delay(:,:,3) = angle(d.normData);
    d2.amp(:,:,1) = abs(d.xData(1:end/2,:));
    d2.amp(:,:,2) = abs(d.xData(end/2+1:end,:));
    d2.amp(:,:,3) = abs(d.normData);
    switch p.keepInfo
        case {'m','p','mp'}
        case {'r','i'}
            %remove voxel mean phase estimated from the plaid
            d2.delay(:,:,1) = wrapToPi(d2.delay(:,:,1)-repmat(circ_mean(d2.delay(:,:,3),[],1),[size(d2.delay,1) 1]));
            d2.delay(:,:,2) = wrapToPi(d2.delay(:,:,2)-repmat(circ_mean(d2.delay(:,:,3),[],1),[size(d2.delay,1) 1]));
            d2.delay(:,:,3) = wrapToPi(d2.delay(:,:,3)-repmat(circ_mean(d2.delay(:,:,3),[],1),[size(d2.delay,1) 1]));
            switch p.keepInfo
                case 'r'
                    %set amplitude to the real value (amplitude at the mean delay)
                    [d2.amp,~] = pol2cart(d2.delay,d2.amp);
                case 'i'
                    %set amplitude to the mag-normalized imaginary value (delay relative to mean delay)
                    [~,d2.amp] = pol2cart(d2.delay,1);
            end
        otherwise
            error('double-check that')
    end
    
    d2.sessionLabel(:,:,1) = d.sessionLabel(1:end/2,:);
    d2.sessionLabel(:,:,2) = d.sessionLabel(end/2+1:end,:);
    d2.sessionLabel(:,:,3) = d.sessionLabel(1:end/2,:);
    d2.label(:,:,1) = d.label(1:end/2,:);
    d2.label(:,:,2) = d.label(end/2+1:end,:);
    d2.label(:,:,3) = ones(size(d.label(1:end/2,:))).*3;
    if ~strcmp(p.k,'randK')
        d2.crossVal(:,:,1) = d.crossVal(1:end/2,:);
        d2.crossVal(:,:,2) = d.crossVal(end/2+1:end,:);
        d2.crossVal(:,:,3) = nan(size(d.crossVal(1:end/2,:)));
    end
    
    switch p.k
        %         case {'auto','autoCmplt'}
        %             error('teInd and trInd badly defined')
        %             teInd = squeeze(d2.crossVal(fold,1,1:2));
        %             d2tr = d2;
        %             allFields = fields(d2tr);
        %             for i = 1:length(allFields)
        %                 d2te.(allFields{i}) = d2tr.(allFields{i})(teInd,:,:);
        %                 d2tr.(allFields{i})(teInd,:,:) = [];
        %             end
        case 'autoRun'
            teInd1 = d2.crossVal(:,1,1)==fold;
            teInd2 = d2.crossVal(:,1,2)==fold;
            allFields = fields(d2);
            for i = 1:length(allFields)
                d2te.(allFields{i})(:,:,1) = d2.(allFields{i})(teInd1,:,1);
                d2te.(allFields{i})(:,:,2) = d2.(allFields{i})(teInd2,:,2);
                tmp1 = d2.(allFields{i})(~teInd1,:,1);
                tmp2 = d2.(allFields{i})(~teInd2,:,2);
                d2tr.(allFields{i}) = cat(3,tmp1,tmp2); clear tmp1 tmp2
            end
            %         case 'sessBal'
            %             error('double-check that')
            %         case 'randK'
            %             % Define train and test data (purely random selection of one sample per session for the testing set, puts more weights on sessions with fewer sample)
            %             sessions = unique(d2.sessionLabel(:,:,1));
            %             for i = 1:length(sessions)
            %                 curInd = find(d2.sessionLabel(:,:,1)==sessions(i));
            %                 teInd(i) = curInd(randperm(length(curInd),1));
            %             end
            %             d2tr = d2;
            %             allFields = fields(d2tr);
            %             for i = 1:length(allFields)
            %                 d2te.(allFields{i}) = d2tr.(allFields{i})(teInd,:,:);
            %                 d2tr.(allFields{i})(teInd,:,:) = [];
            %             end
        otherwise
            error('double-check that')
    end
    
    
    if p.regSession
        error('double-check that')
        %estimate session effect from training set
        [betaHat,X] = estimateSessionEffect(d2tr.amp,d2tr.sessionLabel);
        firsSessInd = size(X,2)-(length(unique(d2tr.sessionLabel))-1)+1;
        sessionEffect = betaHat(:,firsSessInd:end)';
        
        %remove session effect
        for sess = 2:length(unique(d2tr.sessionLabel))
            %from training set
            sessInd = d2tr.sessionLabel(:,:,1)==sess;
            d2tr.amp(sessInd,:,:) = d2tr.amp(sessInd,:,:)-repmat(sessionEffect,[length(find(sessInd)) 1 3]);
            %from testing set
            sessInd = d2te.sessionLabel(:,:,1)==sess;
            d2te.amp(sessInd,:,:) = d2te.amp(sessInd,:,:)-repmat(sessionEffect,[length(find(sessInd)) 1 3]);
        end
    end
    
    
    
    %% Normalize
    allFields = fields(d2tr);
    for i = 1:length(allFields)
        tr.(allFields{i}) = cat(1,d2tr.(allFields{i})(:,:,1),d2tr.(allFields{i})(:,:,2));
    end
    allFields = fields(d2te);
    for i = 1:length(allFields)
        te.(allFields{i}) = cat(1,d2te.(allFields{i})(:,:,1),d2te.(allFields{i})(:,:,2));
    end
    trPreNorm = tr;
    tePreNorm = te;
    % Delay
    %remove mean phase and rotate
    phase_shift = circ_mean(tr.delay,[],1)-p.rotTo;
    tr.delay = wrapToPi(tr.delay-repmat(phase_shift,size(tr.delay,1),1));
    te.delay = wrapToPi(te.delay-repmat(phase_shift,size(te.delay,1),1));
    switch p.keepInfo
        case {'m','r','i'}
            % Amp
            %shift
            amp_shift = mean(tr.amp,1);
            tr.amp = tr.amp-repmat(amp_shift,size(tr.amp,1),1);
            te.amp = te.amp-repmat(amp_shift,size(te.amp,1),1);
            %scale
            amp_scale = std(tr.amp,[],1);
            tr.amp = tr.amp./repmat(amp_scale,size(tr.amp,1),1);
            te.amp = te.amp./repmat(amp_scale,size(te.amp,1),1);
        case 'mp'
            %Amp scale only (equivalent to image scaling only)
            amp_scale = mean(tr.amp,1);
            tr.amp = tr.amp./repmat(amp_scale,[size(tr.amp,1) 1]);
            te.amp = te.amp./repmat(amp_scale,[size(te.amp,1) 1]);
            
            %             %Image Shift and Scale
            %             [trR,trI] = pol2cart(tr.delay,tr.amp);
            %             [teR,teI] = pol2cart(te.delay,te.amp);
            %             %Shift
            %             %real
            %             trR_shift = mean(trR,1);
            %             trR = trR-repmat(trR_shift,size(trR,1),1);
            %             teR = teR-repmat(trR_shift,size(teR,1),1);
            %             %imag
            %             trI_shift = mean(trI,1);
            %             trI = trI-repmat(trI_shift,size(trI,1),1);
            %             teI = teI-repmat(trI_shift,size(teI,1),1);
            %             %Scale
            %             tr_scale = mean([std(trR,[],1); std(trI,[],1)],1);
            %             trR = trR./repmat(tr_scale,size(trR,1),1); trI = trI./repmat(tr_scale,size(trI,1),1);
            %             teR = teR./repmat(tr_scale,size(teR,1),1); teI = teI./repmat(tr_scale,size(teI,1),1);
            %
            %             [tr.delay,tr.amp] = cart2pol(trR,trI);
            %             [te.delay,te.amp] = cart2pol(teR,teI);
            %             clear trR trI tr_shift tr_scale teR teI te_shift te_scale
        case {'p'}
        otherwise
            error('did not code that')
    end
    
    switch p.algorithm
        case 'runSVM'
            %% Sorting for feature slection (previously done with newFeatSorting.m)
            % compute sorting
            if ~exist('sortIndFix','var')
                switch p.keepInfo
                    case {'m','r','i'}
                        p.sort.feat = computeMI(tr.amp,tr.label==1,tr.label==2)';
                    case 'p'
                        for dim = 1:size(tr.delay,2)
                            [pval(dim), table] = circ_wwtest(tr.delay(tr.label==1,dim),tr.delay(tr.label==2,dim),[],[],0);
                        end
                        p.sort.feat = 1-pval;
                    case 'mp'
                        [X,Y] = pol2cart(tr.delay,tr.amp);
                        p.sort.feat = computeMI(complex(X,Y),tr.label==1,tr.label==2)';
                    otherwise
                        error('did not code that')
                end
            else
                p.sort.feat = nan(1,p.nFeatures);
            end
            
            % Sorting for functional ROI
            p.sort.fun = p.funcROI.vec.stats;
            [~,p.sort.indFun] = sort(p.sort.fun,'descend');
            
            trPreFunSort = tr;
            tePreFunSort = te;
            
            
            
            
            %% Loop over functional ROI selection levels
            for funLevInd = 1:length(p.featLevList)
                % Apply fun sorting (and select fun)
                trFeatSort.amp = trPreFunSort.amp(:,p.sort.indFun((1:p.featLevList(funLevInd))));
                trFeatSort.delay = trPreFunSort.delay(:,p.sort.indFun((1:p.featLevList(funLevInd))));
                teFeatSort.amp = tePreFunSort.amp(:,p.sort.indFun((1:p.featLevList(funLevInd))));
                teFeatSort.delay = tePreFunSort.delay(:,p.sort.indFun((1:p.featLevList(funLevInd))));
                curFun = p.sort.fun(p.sort.indFun((1:p.featLevList(funLevInd))));
                curFeat = p.sort.feat(p.sort.indFun((1:p.featLevList(funLevInd))));
                
                % Sort feat
                [~,curIndFeat] = sort(curFeat,'descend');
                
                % Apply feat sorting
                trFeatSort.amp = trFeatSort.amp(:,curIndFeat);
                trFeatSort.delay = trFeatSort.delay(:,curIndFeat);
                teFeatSort.amp = teFeatSort.amp(:,curIndFeat);
                teFeatSort.delay = teFeatSort.delay(:,curIndFeat);
                curFun = curFun(curIndFeat);
                curFeat = curFeat(curIndFeat);
                
                %% Loop over feature selection levels
                for featLevInd = 1:length(p.featLevList)
                    if p.featLevList(featLevInd)>length(curFun)
                        break
                    end
                    % Select feat
                    tr.amp = trFeatSort.amp(:,1:p.featLevList(featLevInd));
                    tr.delay = trFeatSort.delay(:,1:p.featLevList(featLevInd));
                    te.amp = teFeatSort.amp(:,1:p.featLevList(featLevInd));
                    te.delay = teFeatSort.delay(:,1:p.featLevList(featLevInd));
                    curFun(1:p.featLevList(featLevInd));
                    curFeat(1:p.featLevList(featLevInd));
                    
                    % Run SVM
                    %Train and test
                    switch p.keepInfo
                        case {'m','r','i'}
                            complexFlag = 0;
                            trFinal.data = tr.amp;
                            trFinal.label = tr.label;
                            teFinal.data = te.amp;
                            teFinal.label = te.label;
                        case 'p'
                            complexFlag = 1;
                            [X,Y] = pol2cart(tr.delay,1);
                            trFinal.data = [X Y];
                            trFinal.label = tr.label;
                            [X,Y] = pol2cart(te.delay,1);
                            teFinal.data = [X Y];
                            teFinal.label = te.label;
                        case 'mp'
                            complexFlag = 1;
                            [X,Y] = pol2cart(tr.delay,tr.amp);
                            trFinal.data = [X Y];
                            trFinal.label = tr.label;
                            [X,Y] = pol2cart(te.delay,te.amp);
                            teFinal.data = [X Y];
                            teFinal.label = te.label;
                    end
                    
                    [hitRate(fold,funLevInd,featLevInd),~] = trainNtestSVM(trFinal.data,trFinal.label,teFinal.data,teFinal.label,p);
                    
                    %                 %Extract w
                    %                 alpha = zeros(length(svmStruct.GroupNames),1);
                    %                 alpha(svmStruct.SupportVectorIndices,1) = svmStruct.Alpha;
                    %                 y = (svmStruct.GroupNames-1)*2-1;
                    %                 x = trFinal.data;
                    %                 wTmp = zeros(1,size(x,2));
                    %                 for i = 1:size(x,1)
                    %                     wTmp = wTmp + alpha(i)*y(i)*x(i,:);
                    %                 end
                    %
                    %                 %Extract s
                    %                 s = wTmp*x';
                    %
                    %                 %Extract A
                    %                 Atmp = wTmp*cov(x)/cov(s);
                    %
                    %                 %Compile
                    %                 if complexFlag
                    %                     wTmp = complex(wTmp(1:end/2),wTmp(end/2+1:end));
                    %                     Atmp = complex(Atmp(1:end/2),Atmp(end/2+1:end));
                    %                 end
                    %                 wFold(funLevInd,featLevInd,1:p.featLevList(featLevInd)) = wTmp;
                    %                 Afold(funLevInd,featLevInd,1:p.featLevList(featLevInd)) = Atmp;
                end
            end
            %     if p.doSVM
            %         %Sort back weights to original and store
            %         w(:,:,fold) = sort_back(wFold,repmat(sortInd(fold,:),size(wFold,1),1),2);
            %         A(:,:,fold) = sort_back(Afold,repmat(sortInd(fold,:),size(Afold,1),1),2);
            %         %         data(:,:,fold) = sort_back(dataFold,repmat(sortInd(fold,:),size(dataFold,1),1),2);
            %     end
            if verbose
                toc
            end
        case 'runSVM_RFE'
            switch p.keepInfo
                case {'m','r','i'}
                    complexFlag = 0;
                    trFinal.data = tr.amp;
                    trFinal.label = tr.label;
                    teFinal.data = te.amp;
                    teFinal.label = te.label;
                case 'p'
                    complexFlag = 1;
                    [X,Y] = pol2cart(tr.delay,1);
                    trFinal.data = [X Y];
                    trFinal.label = tr.label;
                    [X,Y] = pol2cart(te.delay,1);
                    teFinal.data = [X Y];
                    teFinal.label = te.label;
                case 'mp'
                    complexFlag = 1;
                    [X,Y] = pol2cart(tr.delay,tr.amp);
                    trFinal.data = [X Y];
                    trFinal.label = tr.label;
                    [X,Y] = pol2cart(te.delay,te.amp);
                    teFinal.data = [X Y];
                    teFinal.label = te.label;
            end
            
            tmpFeatLevList = sort(p.featLevList,'descend');
            if tmpFeatLevList(1)~=p.nFeatures
                error('last feat level must contain all features')
            end
            
            %Run RFE R times
            for R = 1:p.RFE_R
                trFinalR = trFinal;
                teFinalR = teFinal;
                ind2keep = 1:p.nFeatures;
                %Define L splits
                LsplitInd = [];
                for L = 1:p.RFE_L
                    LsplitInd = [LsplitInd [randperm(size(trFinalR.data,1)/2)'>round(size(trFinalR.data,1)/2*(1/3)); randperm(size(trFinalR.data,1)/2)'>size(trFinalR.data,1)/2*(1/3)]];
                end
                LsplitInd = logical(LsplitInd);
                for featLevel = 1:length(tmpFeatLevList)
                    %Get w on L splits
                    w = nan(p.RFE_L,p.nFeatures);
                    for L = 1:p.RFE_L
                        %split
                        trFinalL = trFinalR;
                        trFinalL.data = trFinalL.data(LsplitInd(:,L),:);
                        trFinalL.label = trFinalL.label(LsplitInd(:,L),:);
                        
                        %train and assess generalization
%                         [hitRate(fold,R,featLevel,L),svmStruct] = trainNtestSVM(trFinalL.data,trFinalL.label,teFinalR.data,teFinalR.label,p);
                        [~,svmStruct] = trainNtestSVM(trFinalL.data,trFinalL.label,teFinalR.data,teFinalR.label,p);
                        %w
                        alpha = zeros(size(trFinalL.data,1),1);
                        alpha(svmStruct.SupportVectorIndices,1) = svmStruct.Alpha;
                        y = (trFinalL.label-1)*2-1;
                        x = trFinalL.data;
                        curw = zeros(1,size(x,2));
                        for i = 1:size(x,1)
                            curw = curw + alpha(i)*y(i)*x(i,:);
                        end
                        if complexFlag
                            [~,w(L,ind2keep)] = cart2pol(curw(1:end/2),curw(end/2+1:end));
%                             [~,w(fold,R,featLevel,L,ind2keep)] = cart2pol(curw(1:end/2),curw(end/2+1:end));
                        else
                            w(L,ind2keep) = curw;
%                             w(fold,R,featLevel,L,ind2keep) = curw;
                        end
                    end
                    %Generalization performance
                    [hitRate(fold,R,featLevel),~] = trainNtestSVM(trFinalR.data,trFinalR.label,teFinalR.data,teFinalR.label,p);
                    %Eliminate features
                    if featLevel~=length(tmpFeatLevList)
                        validInd = find(~isnan(w(1,:)));
                        curw = squeeze(w(:,validInd));
%                         validInd = find(~isnan(w(fold,R,featLevel,1,:)));
%                         curw = squeeze(w(fold,R,featLevel,:,validInd));
                        [~,b] = sort(mean(abs(curw),1),'descend');
                        ind2keep = sort(validInd(b(1:tmpFeatLevList(featLevel+1))));
                        if complexFlag
                            ind2keep = [ind2keep ind2keep+size(trFinal.data,2)/2];
                        end
                        trFinalR.data = trFinal.data(:,ind2keep);
                        teFinalR.data = teFinal.data(:,ind2keep);
                        if complexFlag
                            ind2keep = ind2keep(1:end/2);
                        end
                    end
                end
            end
        otherwise
            error('did not code that')
    end
end

%Average over cross-validation folds
switch p.algorithm
    case 'runSVM'
        hitRate = mean(hitRate,1);
    case 'runSVM_RFE'
%         hitRate = mean(mean(mean(hitRate,1),2),4);
        hitRate = mean(mean(hitRate,1),2);
end


