function [resBS,resBShr,resWS,f] = runDecoding(p,verbose,paths,resPall,featSelPall)
if ~exist('verbose','var')
    verbose = 1;
end
% if ~exist('figOption','var') || isempty(figOption)
%     figOption.save = 0;
%     figOption.subj = 1; % 'all' or subjInd
% end
if ~isfield(p,'chanSpace') || isempty(p.chanSpace)
    p.chanSpace = 'cart'; % 'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'
    % 'hr' 'hrNoAmp' 'cart' 'cartNoAmp' cartNoAmp_HT 'cartReal', 'cartImag', 'pol', 'polMag' 'polMag_T' or 'polDelay'
end
if ~exist('resPall','var') || isempty(resPall)
    permFlag = 0;
    resPall = [];
    featSelPall = [];
else
    permFlag = 1;
end
f = [];
% tic

%% Define paths
subjList = paths.subjList;
% repoPath = paths.repoPath;
    funPath = paths.funPath;
        inDir  = paths.inDir;
        inDir2  = paths.inDir2;
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'funPath' 'inDir' 'inDir2'}
% for tmp = {'repoPath' 'funPath' 'inDir' 'inDir2'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp



%% Load response data
dAll = cell(size(subjList,1),1);
for subjInd = 1:length(subjList)
    curFile = fullfile(funPath,'resp',[subjList{subjInd} '.mat']);
    if verbose; disp(['loading: ' curFile]); end
    if ~permFlag
        load(curFile,'resp');
    else
        resp = resPall{subjInd}; resPall{subjInd} = {};
%         load(curFile,'resP');
%         res = resP; clear resP
    end
    dAll{subjInd} = resp; clear resp
end
d = dAll; clear dAll
sessList = fields(d{1});
% Load feature slection
if ~permFlag
    load(fullfile(funPath,'featSel.mat'),'featSel');
else
    featSel = featSelPall; clear featSelPall
%     load(fullfile(funPath,inDir2,'featSel.mat'),'featSelP');
%     featSel = featSelP; clear featSelP
end
if verbose
    disp('---');
    disp(['Channel space: ' p.chanSpace '-' p.condPair]);
    disp(['Complex space: ' p.complexSpace]);
end

%% Reorganize
dP = cell(size(d,2),length(sessList));
for subjInd = 1:length(d)
    for sessInd = 1:length(sessList)
        dP{subjInd,sessInd} = d{subjInd}.(sessList{sessInd});
        d{subjInd}.(sessList{sessInd}) = [];
        dP{subjInd,sessInd}.featSel = featSel{subjInd,sessInd};
    end
end
d = dP; clear dP

% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc

%% Define feature selection
if ~strcmp(featSel{1,1}.featSeq.info2,p.featSel.global.method)
    error(['p.featSel.global.method and featSel.featSeq.info2 not matching' newline 'Try reruning processFeatSel.m'])
end
featSelSteps_labelList = featSel{1,1}.featSeq.featSelList;
featSelConds_labelList = featSel{1,1}.featSeq.condPairList;
[ind_nSpecFeatSel,ind_nSpecFeatSelCond,ind_specFeatSel,ind_specFeatSelCond] = defineFeatSel(featSelSteps_labelList,featSelConds_labelList,p.featSel.global.method,p.condPair);
for subjInd = 1:size(d,1)
    for sessInd = 1:size(d,2)
        featSel{subjInd,sessInd}.indIn = ...
            all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond),2)...
            & all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_specFeatSel,ind_specFeatSelCond),2);
    end
end

% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc

%% Example plot of trigonometric (polar) representation
subjInd = p.figOption.subjInd;
sessInd = p.figOption.sessInd;
if p.figOption.verbose==1
    f = [f plotNorm(d{subjInd,sessInd},p,featSel{subjInd,sessInd},[],0)];
elseif p.figOption.verbose>1
    f = [f plotNorm(d{subjInd,sessInd},p,featSel{subjInd,sessInd},[],1)];
end
% if p.figOption.verbose>=1 && figOption.save
%     error('code that')
%     filename = fullfile(pwd,mfilename);
%     if ~exist(filename,'dir'); mkdir(filename); end
%     filename = fullfile(filename,p.chanSpace);
%     f.Color = 'none';
%     set(findobj(f.Children,'type','Axes'),'color','none')
%     set(findobj(f.Children,'type','PolarAxes'),'color','none')
%     saveas(f,[filename '.svg']); if verbose; disp([filename '.svg']); end
%     f.Color = 'w';
%     set(findobj(f.Children,'type','Axes'),'color','w')
%     set(findobj(f.Children,'type','PolarAxes'),'color','w')
%     saveas(f,[filename '.fig']); if verbose; disp([filename '.fig']); end
%     saveas(f,[filename '.jpg']); if verbose; disp([filename '.jpg']); end
% end

% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc


if p.svm.doWithin
    %% Within-session SVM cross-validation (with cross-session feature selection)
    svmModelK = cell(size(d));
    nrmlzK = cell(size(d));
    yHatK = cell(size(d));
    yHatK_tr = cell(size(d));
    Y = cell(size(d));
    K = cell(size(d));
    Yhr = cell(size(d));
    Khr = cell(size(d));
    nOrig = nan(size(d));
    n = nan(size(d));
    
    if verbose
        disp('Within-session SVM cross-validation')
        disp('with between-session feature selection')
    end
    for i = 1:numel(d)
        if verbose
            disp(['-training and testing sess ' num2str(i) '/' num2str(numel(d))])
        end
        [subjInd,sessInd] = ind2sub(size(d),i);
        switch sessInd
            case 1
                sessIndCross = 2;
            case 2
                sessIndCross = 1;
            otherwise
                error('X')
        end
        
        % get data
        [X,y,k] = getXYK(d{subjInd,sessInd},p);
        
        % cross-session feature selection
        featSelInd = featSel{subjInd,sessIndCross}.indIn;
        x = X(:,featSelInd);
        
        % k-fold training
        [svmModelK{i},nrmlzK{i}] = SVMtrain(y,x,p,k);
        % k-fold testing
        [yHatK{i},yHatK_tr{i}] = SVMtest(y,x,svmModelK{i},nrmlzK{i},k);
        yHatK{i} = nanmean(yHatK{i},2);
        % output
        Y{i} = y;
        K{i} = k;
        nOrig(i) = length(featSelInd);
        n(i) = nnz(featSelInd);
    end

    
    % line_num=dbstack; % get this line number
    % disp(['line ' num2str(line_num(1).line)]) % displays the line number
    % toc
end

%% Within-session SVM training (and feature selection) + Cross-session SVM testing
if verbose
    disp('Between-session SVM cross-validation')
    disp('with feature selection from the training session')
end
svmModel = cell(size(d));
nrmlz = cell(size(d));
yHat_tr = cell(size(d));
yHat = cell(size(d));
yHatHr = cell(size(d));
yHatHr_nSpec = cell(size(d));
if ~p.svm.doWithin
    Y = cell(size(d));
    K = cell(size(d));
    Yhr = cell(size(d));
    Khr = cell(size(d));
    nOrig = nan(size(d));
    n = nan(size(d));
end
% Training
for i = 1:numel(d)
    if verbose
        disp(['-training sess ' num2str(i) '/' num2str(numel(d))])
    end
    [subjInd,trainInd] = ind2sub(size(d),i);
    % get training data
    [X,y,k] = getXYK(d{subjInd,trainInd},p);
    % same-session feature-selection
    featSelInd = featSel{subjInd,trainInd}.indIn;
    x = X(:,featSelInd);
    % training
    [svmModel{subjInd,trainInd},nrmlz{subjInd,trainInd}] = SVMtrain(y,x,p);
    % same session testing (not crossvalidated)
    [yHat_tr{subjInd,trainInd},~] = SVMtest(y,x,svmModel{subjInd,trainInd},nrmlz{subjInd,trainInd});
    
    Y{i} = y;
    K{i} = k;
    nOrig(i) = length(featSelInd);
    n(i) = nnz(featSelInd);
end

% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc

% Testing
for i = 1:numel(d)
    if verbose
        disp(['-testing sess ' num2str(i) '/' num2str(numel(d))])
    end
    [subjInd,testInd] = ind2sub(size(d),i);
    switch testInd
        case 1
            trainInd = 2;
        case 2
            trainInd = 1;
        otherwise
            error('X')
    end
    % get test data
    [X,y,~] = getXYK(d{subjInd,testInd},p);
    % cross-session feature selection
    featSelInd = featSel{subjInd,trainInd}.indIn;
    x = X(:,featSelInd);
    % cross-session SVM testing
    [yHat{subjInd,testInd},~] = SVMtest(y,x,svmModel{subjInd,trainInd},[],[]);
end

% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc

if ~permFlag...
        && strcmp(p.svm.kernel.type,'lin')...
        && strcmp(p.svm.complexSpace,'bouboulisDeg1')
    %% Extract channel HR
    % Testing
    for i = 1:numel(d)
        if verbose
            disp(['-channel timeseries for sess ' num2str(i) '/' num2str(numel(d))])
        end
        [subjInd,testInd] = ind2sub(size(d),i);
        switch testInd
            case 1
                trainInd = 2;
            case 2
                trainInd = 1;
            otherwise
                error('X')
        end
        % get test data
        pTmp = p;
        pTmp.condPair = 'all';
        hrFlag = true;
        [x,y,k] = getXYK(d{subjInd,testInd},pTmp,hrFlag);
        % cross-session feature selection
        featSelInd = featSel{subjInd,trainInd}.indIn;
        x = x(:,featSelInd,:);
        % transform w
        curSvmModel = svmModel{subjInd,trainInd};
        switch p.chanSpace
            case 'cart'
                % complexify x?
                % realifiy w?
                
                % Quadrant-specific analysis:
                % 1) +real
                % 2) +real +imag
                % 3) +imag
                % 4) +imag -real
                w = complex(curSvmModel.w(1:end/2),curSvmModel.w(end/2+1:end));
                %             figure('WindowStyle','docked');
                %             polarplot(w,'k.'); hold on
                secAngle = 0:pi/4:pi/4*3;
                wSec = zeros(length(secAngle),length(w));
                for secInd = 1:length(secAngle)
                    indA = (wrapToPi(angle(w)-secAngle(secInd)) < pi/4/2 ...
                        & wrapToPi(angle(w)-secAngle(secInd)) >= -pi/4/2);
                    indB = (wrapToPi(angle(w)-secAngle(secInd) - pi) < pi/4/2 ...
                        & wrapToPi(angle(w)-secAngle(secInd) - pi) >= -pi/4/2);
                    %                 polarplot(w(indA | indB),'.'); hold on
                    wSec(secInd,indA) = abs(w(1,indA));
                    wSec(secInd,indB) = -abs(w(1,indB));
                end
                curSvmModel.w = wSec;
            case 'cartNoAmp'
                % complexify w
                w = complex(curSvmModel.w(1:end/2),curSvmModel.w(end/2+1:end));
                % defining the sign of wi according to its angle is what allows
                % it to oppose the conditions to be compared a wi complex
                % vector becomes its abolute value of sign determined by tha
                % angle of the vector
                curSvmModel.w = sign(angle(w)).*abs(w);
            case 'cartNoDelay'
                % nothing to do here, everything is real valued (originating in the abs value of xi reponse vectors)
            case 'delay'
                % nothing to do here, everything is real valued (originating in the angle value of xi reponse vectors)
                % BUT note that this angle is stored in the imaginary part
                % of the xi
            otherwise
                error('X')
        end
        %     if length(curSvmModel.w)==2*size(x,2)
        %         curSvmModel.w = abs(complex(curSvmModel.w(1:end/2),curSvmModel.w(end/2+1:end)));
        %     end
        % get channel response
        for tInd = 1:size(x,3)
            % cross-session SVM testing
            [yHatHr{subjInd,testInd}(:,:,tInd,:),~] = SVMtest(y,x(:,:,tInd),curSvmModel,0,[]);
        end
        % get non-specific channel response
        curSvmModel.w = abs(curSvmModel.w);
        for tInd = 1:size(x,3)
            % cross-session SVM testing
            [yHatHr_nSpec{subjInd,testInd}(:,:,tInd,:),~] = SVMtest(y,x(:,:,tInd),curSvmModel,0,[]);
        end
        % output
        Yhr{i} = y;
        Khr{i} = k;
    end
end


%% Compute performance metrics (SLOW 10sec/16sec)
% Between-sess
if verbose
    disp('BS perfMetric: computing')
end
resBS_sess = repmat(perfMetric,[size(d,1) size(d,2)]);
if permFlag
    for i = 1:size(d,1)
        k = cat(1,K{i,1},K{i,2}+max(K{i,1}));
        y = cat(1,Y{i,:});
        yhat = cat(1,yHat{i,:});
        resBS_sess(i) = perfMetric(y,yhat,k,1);
    end
    resBS.auc = nan(size(resBS_sess,1),1);
    resBS.auc(:) = [resBS_sess.auc];
    resBShr = [];
    resWS = [];
    f = [];
    return
else
    for i = 1:numel(d)
        resBS_sess(i) = perfMetric(Y{i},yHat{i},K{i});
    end
end
acc_fdr = mafdr([resBS_sess.acc_p],'BHFDR',true);
for i = 1:numel(d)
    resBS_sess(i).acc_fdr = acc_fdr(i);
    resBS_sess(i).nVoxOrig = size(d{i}.sin,1);
    resBS_sess(i).nVox = nnz(featSel{i}.indIn);
    resBS_sess(i).chanSpace = p.chanSpace;
    resBS_sess(i).complexSpace = p.complexSpace;
    resBS_sess(i).svmKernel = p.svm.kernel.type;
    resBS_sess(i).condPair = p.condPair;
end
% resBS_sess = orderfields(resBS_sess,[1 2 3 4 5 6 7 8 9 15 10 11 12 13 14 16 17 18 19]);
% resBS_sess = orderfields(resBS_sess,[1 2 3 4 5 6 7 8 9 17 10 11 12 13 14 15 16 18 19 20 21]);
disp('BS perfMetric: done')

if ~isempty(Yhr{1})
% Between-sess Channel HR
resBShr_sess = repmat(perfMetric,[size(d,1) size(d,2)]);
disp('BShr perfMetric: computing')
for i = 1:numel(d)
    disp(['sess' num2str(i) '/' num2str(numel(d))])
    resBShr_sess(i) = perfMetric(Yhr{i},yHatHr{i},Khr{i});
end
for i = 1:numel(d)
    resBShr_sess(i).yHat_nSpec = {permute(yHatHr_nSpec{i},[1 3 4 2])};
end
for i = 1:numel(d)
    resBShr_sess(i).acc_fdr = [];
    resBShr_sess(i).nVoxOrig = size(d{i}.sin,1);
    resBShr_sess(i).nVox = nnz(featSel{i}.indIn);
    resBShr_sess(i).chanSpace = p.chanSpace;
    resBShr_sess(i).condPair = p.condPair;
end
% resBShr_sess = orderfields(resBShr_sess,[1 2 3 4 5 6 7 8 9 17 10 11 12 13 14 15 16 18 19 20 21 22]);
disp('BShr perfMetric: done')
end

if p.svm.doWithin
    % Within-sess
    disp('WS perfMetric: computing')
    resWS_sess = repmat(perfMetric,[size(d,1) size(d,2)]);
    for i = 1:numel(d)
        resWS_sess(i) = perfMetric(Y{i},yHatK{i},K{i});
    end
    acc_fdr = mafdr([resWS_sess.acc_p],'BHFDR',true);
    for i = 1:numel(d)
        resWS_sess(i).acc_fdr = acc_fdr(i);
        resWS_sess(i).nVoxOrig = size(d{i}.sin,1);
        resWS_sess(i).nVox = nnz(featSel{i}.indIn);
        resWS_sess(i).chanSpace = p.chanSpace;
        resWS_sess(i).complexSpace = p.complexSpace;
        resWS_sess(i).svmKernel = p.svm.kernel.type;
        resWS_sess(i).condPair = p.condPair;
    end
    % resWS_sess = orderfields(resWS_sess,[1 2 3 4 5 6 7 8 9 17 10 11 12 13 14 15 16 18 19 20 21]);
    disp('WS perfMetric: done')
end

% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc

%% Summarize group performances (SLOW 4sec/16sec)
disp('perfMetric summary: computing')
[resBSsess,resBSsubj,resBSgroup] = summarizePerf(resBS_sess);
if exist('resBShr_sess','var')
    [resBShrSess,~,~] = summarizePerf(resBShr_sess);
end

if p.svm.doWithin
    [resWSsess,resWSsubj,resWSgroup] = summarizePerf(resWS_sess);
    disp('perfMetric summary: computing')
end

% % Combine WS and BS
% fieldList = fields(resBSsubj);
% for i = 1:length(fieldList)
%     if isnumeric(resBSsubj.(fieldList{i}))
%         resSubj.(fieldList{i}) = mean(cat(3,resWSsubj.(fieldList{i}),resBSsubj.(fieldList{i})),3);
%     end
% end
% [~,P,~,STATS] = ttest(resSubj.acc,0.5,'tail','right');
% resGroup.acc_T = STATS.tstat;
% resGroup.acc_P = P;
% [~,P,~,STATS] = ttest(resSubj.auc,0.5,'tail','right');
% resGroup.auc_T = STATS.tstat;
% resGroup.auc_P = P;
% [P,~,STATS] = signrank(resSubj.acc,0.5,'tail','right');
% resGroup.acc_wilcoxonSignedrank = STATS.signedrank;
% resGroup.acc_wilcoxonP = P;
% [P,~,STATS] = signrank(resSubj.auc,0.5,'tail','right');
% resGroup.auc_wilcoxonSignedrank = STATS.signedrank;
% resGroup.auc_wilcoxonP = P;


% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc

%% Print info and output
if verbose
    disp('-----------------')
    if p.svm.doWithin
        disp('*Within-session*')
        printRes2(resWSgroup)
        disp(' ')
        %     disp('*Within+Between*')
        %     printRes2(resGroup)
        %     disp(' ')
        disp('*Between-session*')
        printRes2(resBSgroup)
    end
    disp('-----------------')
end

resBS.sess = resBSsess;
resBS.sess.subjList = subjList;
resBS.subj = resBSsubj;
resBS.subj.subjList = subjList;
resBS.group = resBSgroup;

if exist('resBShrSess','var')
    resBShr.sess = resBShrSess;
    resBShr.sess.subjList = subjList;
else
    resBShr = [];
end

if p.svm.doWithin
    resWS.sess = resWSsess;
    resWS.sess.subjList = subjList;
    resWS.subj = resWSsubj;
    resWS.subj.subjList = subjList;
    resWS.group = resWSgroup;
else
    resWS = [];
end

% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc

if p.svm.doWithin
    %% Plot between-session vs within-session
    if verbose
        % Decoding
        sz = size(resWS.sess.acc);
        resBS_tmp = resBS;
        fieldList1 = fields(resBS_tmp);
        for i1 = 1:length(fieldList1)
            fieldList2 = fields(resBS_tmp.(fieldList1{i1}));
            for i2 = 1:length(fieldList2)
                tmp = resBS_tmp.(fieldList1{i1}).(fieldList2{i2});
                if all(size(tmp) == sz)
                    resBS_tmp.(fieldList1{i1}).(fieldList2{i2}) = tmp(:,[2 1]);
                end
            end
        end
        
        if p.figOption.verbose>1
            f = [f figure('WindowStyle','docked','visible','on')];
        else
            f = [f figure('WindowStyle','docked','visible','off')];
        end
        compareRes(resBS,resWS)
        ax = gca;
        ax.XLabel.String = {'between-session' ax.XLabel.String};
        ax.YLabel.String = {'within-session' ax.YLabel.String};
        textLine = {};
        textLine = [textLine {[p.chanSpace '; ' p.condPair]}];
        if p.figOption.verbose>1
            textLine = [textLine {strjoin(featSel{p.figOption.subjInd,1}.featSeq.featSelList,'; ')}];
        end
        
        if p.figOption.verbose>2
            n    = nan([length(featSel{1}.featSeq.featSelList)+1 size(featSel)]);
            %     n = nan([length(featSel{1}.n) size(featSel)]);
            for sessInd = 1:numel(featSel)
                tmp = nan(size(featSel{sessInd}.featSeq.featIndIn(:,:,1)));
                tmp(:,ind_nSpecFeatSel) = featSel{sessInd}.featSeq.featIndIn(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond);
                tmp(:,ind_specFeatSel) = featSel{sessInd}.featSeq.featIndIn(:,ind_specFeatSel,ind_specFeatSelCond);
                n(2:end,sessInd) = sum(cumsum(tmp,2)==1:length(ind_nSpecFeatSel),1);
                %         n(2:end,sessInd) = sum(cumsum(featSel{sessInd}.featSeq.featIndIn(:,:,featSel{sessInd}.condPairCurInd),2) == 1:size(featSel{sessInd}.featSeq.featIndIn,2),1);
                n(1,sessInd) = size(tmp,1);
            end
            nMin = min(n(:,:),[],2)';
            nMin = cellstr(num2str(nMin'))';
            nMean = round(mean(n(:,:),2)');
            nMean = cellstr(num2str(nMean'))';
            nMax = max(n(:,:),[],2)';
            nMax = cellstr(num2str(nMax'))';
            
            textLine = [textLine {['min: ' strjoin(nMin,'; ')]}];
            textLine = [textLine {['mean: ' strjoin(nMean,'; ')]}];
            textLine = [textLine {['max: ' strjoin(nMax,'; ')]}];
        end
        
        textLine = strjoin(textLine,newline);
        title(textLine);
        uistack(patch([0 0 0.5 0.5 0],[0.5 0 0 0.5 0.5],[1 1 1].*0.7),'bottom')
    end
end


% line_num=dbstack; % get this line number
% disp(['line ' num2str(line_num(1).line)]) % displays the line number
% toc


function [ind2,info] = getFeatInd(featMap,featSelInfo,method,ind,dir)
switch method
    case '%ile'
        thresh = prctile(featMap(ind),featSelInfo.(method));
        info = ['prctile' dir num2str(featSelInfo.(method))];
    case 'p'
    case 'fdr'
end
switch dir
    case '<'
        ind2 = featMap<=thresh;
    case '>'
        ind2 = featMap>=thresh;
end
info = [method dir num2str(featSelInfo.(method))];




function [x,y,k,t] = getXYK_wave(dP,p,opt)
% Define x(data), y(label) and k(xValFolds)
switch p.chanSpace
    case {'cart' 'cart_HT' 'cart_HTbSess'...
            'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'...
            'cartNoAmpImag'...
            'cartNoDelay' 'cartNoDelay_HT' 'cartNoDelay_HTbSess'...
            'cartReal' 'cartReal_T'...
            'cartImag' 'cartReal_T'...
            'polMag' 'polMag_T'...
            'polDelay'}
        dP.wave(:,:,:,~dP.good) = [];
        ptsPerCycle = 12;
        
        switch opt
            case 'full'
                x = permute(dP.wave,[2 3 4 1]);
            case 'trial'
                sz = size(dP.wave);
                sz = [sz(1:3) ptsPerCycle sz(4)/ptsPerCycle];
                x = permute(mean(reshape(dP.wave,sz),4),[2 3 5 1 4]);
            case 'trialSparse'
                x = permute(dP.wave(:,:,:,1:ptsPerCycle:end),[2 3 4 1]);
            case 'trialSparseCat2'
                sz = size(dP.wave);
                sz(4) = sz(4)./ptsPerCycle;
                x = nan([sz([2 3 4 1]) ptsPerCycle]);
                for t = 1:ptsPerCycle
                    x(:,:,:,:,t) = permute(dP.wave(:,:,:,t:ptsPerCycle:end),[2 3 4 1]);
                end
                x = permute(x,[2 3 4 1 5]);
                x = permute(x(:,:,:,:),[4 1 2 3]);
            case 'run'
                x = mean(dP.wave,4);
                x = permute(x,[2 3 4 1]);
            otherwise
                error('X')
        end
        sz = size(x);
        y = repmat(1:sz(2),[1 1 sz(3:4)]);
        t = repmat(permute(1:sz(3),[1 3 2 4]),[1 sz(2) 1 sz(4)]);
        k = repmat(permute(1:sz(4),[1 3 4 2]),[1 sz(2:3) 1]);
        x = x(:,:,:);
        y = y(:,:,:);
        t = t(:,:,:);
        k = k(:,:,:);
        x = cat(1,permute(x(:,1,:),[3 1 2]),permute(x(:,2,:),[3 1 2]));
        y = cat(1,permute(y(:,1,:),[3 1 2]),permute(y(:,2,:),[3 1 2]));
        t = cat(1,permute(t(:,1,:),[3 1 2]),permute(t(:,2,:),[3 1 2]));
        k = cat(1,permute(k(:,1,:),[3 1 2]),permute(k(:,2,:),[3 1 2]));
    otherwise
        error('X')
end


function [svmModel,normTr] = SVMtrain(Y,X,p,K,T)
if exist('T','var')
    tList = unique(T);
    for i = 1:length(tList)
        ind = tList(i)==T;
        x = X(ind,:);
        y = Y(ind,:);
        if isempty(K)
            k = [];
        else
            k = K(ind,:);
        end
        [svmModel(:,:,:,i),normTr(:,:,:,i)] = SVMtrain(y,x,p,k);
    end
    for i = 1:size(svmModel,3)
        svmModel(:,:,i,1).w = mean(cat(1,svmModel(:,:,i,:).w),1);
        svmModel(:,:,i,1).b = mean(cat(1,svmModel(:,:,i,:).b),1);
        tmp = cat(1,normTr(:,:,i,:).pol);
        normTr(:,:,i,1).pol.rhoScale = mean(cat(1,tmp.rhoScale),1);
        normTr(:,:,i,1).pol.thetaShift = mean(cat(1,tmp.thetaShift),1);
        tmp = cat(1,normTr(:,:,i,:).svm);
        normTr(:,:,i,1).svm.scale = mean(cat(1,tmp.scale),1);
        normTr(:,:,i,1).svm.shift = mean(cat(1,tmp.shift),1);
    end
    svmModel = svmModel(:,:,:,1);
    normTr = normTr(:,:,:,1);
    return
end
if ~exist('K','var') || isempty(K) || all(K(:)==0)
    K = ones(size(Y)); % no crossvalidation
    kList = 0;
else
    kList = unique(K);
end

normTr.k = [];
normTr.chanSpace = '';
normTr.pol = struct('rhoScale',[],'thetaShift',[]);
normTr.svm = struct('scale',[],'shift',[]);
normTr = repmat(normTr,[1 1 length(kList)]);
svmModel = repmat(struct('k',[],'chanSpace',[],'norm',[],'w',[],'b',[],'Paramters',[],'info',[]),[1 1 length(kList)]);
for kInd = 1:length(kList)
    % Split train and test
    te = K==kList(kInd);
    y = Y(~te,:);
    x = X(~te,:);
    
    % Normalize
    [x,normTr(:,:,kInd).pol] = polarSpaceNormalization(x,p.chanSpace);
    [x,normTr(:,:,kInd).svm] = cartSpaceNormalization(x,p.chanSpace);
    [x,~,~] = complex2svm(x,p.chanSpace,p.complexSpace);
    normTr(:,:,kInd).k = kList(kInd);
    normTr(:,:,kInd).chanSpace = p.chanSpace;
    
    % Train
    switch p.svm.kernel.type
        case 'lin'
            model = svmtrain(y,x,'-t 0 -q');
            w = model.sv_coef'*model.SVs;
            b = model.rho;
            toUnitNorm = norm(w);
            w = w./toUnitNorm;
            b = b./toUnitNorm;
            svmModel(:,:,kInd).model = [];
            svmModel(:,:,kInd).info = '   yHat = x*w''-b   ';
        case 'rbf'
            model = svmtrain(y,x,'-t 2 -q');
            w = nan(1,size(x,2));
            b = nan;
            svmModel(:,:,kInd).model = model;
            svmModel(:,:,kInd).info = '   [~,~,yHat] = svmpredict(y,x,model,''-q'')   ';
        otherwise
            error('X')
    end
    % Store
    svmModel(:,:,kInd).k = kList(kInd);
    svmModel(:,:,kInd).w = w;
    svmModel(:,:,kInd).b = b;
    svmModel(:,:,kInd).kernel = p.svm.kernel.type;
    svmModel(:,:,kInd).Paramters = model.Parameters;
    svmModel(:,:,kInd).chanSpace = p.chanSpace;
    svmModel(:,:,kInd).complexSpace = p.complexSpace;
end

function [yHat,yHatTr] = SVMtest(Y,X,svmModel,nrmlz,K)
if ~exist('K','var') || isempty(K) || all(K(:)==0)
    K = ones(size(Y)); % no crossvalidation
end
kList = [svmModel.k];

yHat = nan([size(Y,1) length(kList) size(svmModel(1).w,1)]);
yHatTr = nan([size(Y,1) length(kList) size(svmModel(1).w,1)]);
for kInd = 1:length(kList)
    x = X;
    chanSpace    = svmModel(kInd).chanSpace;
    complexSpace = svmModel(kInd).complexSpace;
    kernel = svmModel(kInd).kernel;
    % Split train and test
    if kList(kInd)==0
        te = true(size(Y));
    else
        te = K==kList(kInd);
    end
    
    % Normalize
    if ~exist('nrmlz','var') || isempty(nrmlz)
        [x,~] = polarSpaceNormalization(x,svmModel(kInd).chanSpace);
        [x,~] = cartSpaceNormalization(x,svmModel(kInd).chanSpace);
    elseif isstruct(nrmlz)
        [x,~] = polarSpaceNormalization(x,nrmlz(kInd),te);
        [x,~] = cartSpaceNormalization(x,nrmlz(kInd),te);
    elseif nrmlz==0
    elseif nrmlz==1
        error('nope')
    else
        error('don''t know what this is')
    end
    [x,~,~] = complex2svm(x,chanSpace,complexSpace);
    
    switch kernel
        case 'lin'
            w = svmModel(kInd).w;
            b = svmModel(kInd).b;
            yHat(te,kInd,:) = x(te,:)*w'-b;
            yHatTr(~te,kInd,:) = x(~te,:)*w'-b;
        case 'rbf'
            svmModel(kInd).info
            [~,~,yHat(te,kInd,:)] = svmpredict(Y(te),x(te,:),svmModel(kInd).model,'-q');
            [~,~,yHatTr(~te,kInd,:)] = svmpredict(Y(~te),x(~te,:),svmModel(kInd).model,'-q');
        otherwise
            error('X')
    end
end


function [x,nVox,nDim] = complex2svm(x,chanSpace,complexSpace)
% Output SVM ready data
if ~isreal(x)
    nVox = size(x,2);
    switch chanSpace
        case {'cart' 'cartNoAmp' 'cartNoAmp_affineRot' 'cartNoAmp_affineRot_affineCart' 'cart_roi' 'cart_affineRot'}
            switch complexSpace
                case 'bouboulisDeg1'
                    x = cat(2,real(x),imag(x));
                case 'bouboulisDeg2'
                    xExp = nan(size(x,1),size(x,2),size(x,2));
                    for runInd = 1:size(x,1)
                        xExp(runInd,:,:) = real(x(runInd,:))' * imag(x(runInd,:));
                    end
                    x = cat(2,real(x),imag(x),xExp(:,:)); clear xExp
                otherwise
                    error('X')
            end
        case {'cartNoDelay' 'cartReal' 'cartReal_affineRot'}
            x = real(x);
        case {'delay' 'cartNoAmpImag' 'cartImag' 'cartImag_affineRot' 'cartNoAmpImag_affineRot'}
            x = imag(x);
        otherwise
            error('X')
    end
    nDim = size(x,2);
else
    nVox = size(x,2);
    nDim = size(x,2);
end


function f = plotPolNormExampleVox(x,SVMspace,b)
f = figure('WindowStyle','docked');

%% Polar Normalization
% Plot before
subplot(2,2,1); clear hPP
% polarplot(angle(mean(x,1)),abs(mean(x,1)),'.k'); hold on
polarplot(angle(x(:)),abs(x(:)),'.k'); hold on
polarplot(angle(x(:,b)),abs(x(:,b)),'.r'); hold on
drawnow
% Normalize
x = polarSpaceNormalization(x,SVMspace);
% Plot after
subplot(2,2,2);
% polarplot(angle(mean(x,1)),abs(mean(x,1)),'.k'); hold on
polarplot(angle(x(:)),abs(x(:)),'.k'); hold on
polarplot(angle(x(:,b)),abs(x(:,b)),'.r'); hold on
drawnow

%% Cartesian Normalization
% Plot before
subplot(2,2,3);
% hScat = scatter(real(mean(x,1)),imag(mean(x,1)),'ko','filled'); hold on
hScat1 = scatter(real(x(:)),imag(x(:)),'ko','filled'); hold on
hScat1.MarkerEdgeColor = 'none';
hScat1.SizeData = hScat1.SizeData/2;
hScat2 = scatter(real(x(:,b)),imag(x(:,b)),'ro','filled'); hold on
hScat2.MarkerEdgeColor = 'none';
hScat2.SizeData = hScat2.SizeData/2;
% alpha(hScat1,0.03)
% alpha(hScat2,0.03)
ax = gca;
ax.DataAspectRatio = [1 1 1];
ax.PlotBoxAspectRatio = [1 1 1];
xLim = xlim;
delta = abs(diff(xLim)).*0.1;
if ~(xLim(1)<0)
    xLim(1) = -delta;
end
if ~(xLim(2)>0)
    xLim(2) = +delta;
end
xlim(xLim)

yLim = ylim;
if ~(yLim(1)<0)
    yLim(1) = -delta;
end
if ~(yLim(2)>0)
    yLim(2) = +delta;
end
ylim(yLim)

uistack(plot([0 0],ylim,'-k'),'bottom');
uistack(plot(xlim,[0 0],'-k'),'bottom');
grid on
title('before cartNorm')
xlabel('real')
ylabel('imag')
drawnow

% Normalize
x = cartSpaceNormalization(x,SVMspace);

%Plot after
subplot(2,2,4);
% hScat = scatter(real(mean(x,1)),imag(mean(x,1)),'ko','filled'); hold on
hScat1 = scatter(real(x(:)),imag(x(:)),'ko','filled'); hold on
hScat1.MarkerEdgeColor = 'none';
hScat1.SizeData = hScat1.SizeData./2;
hScat2 = scatter(real(x(:,b)),imag(x(:,b)),'ro','filled'); hold on
hScat2.MarkerEdgeColor = 'none';
hScat2.SizeData = hScat2.SizeData./2;
% alpha(hScat1,0.03)
ax = gca;
ax.DataAspectRatio = [1 1 1];
ax.PlotBoxAspectRatio = [1 1 1];
xLim = xlim;
delta = abs(diff(xLim)).*0.1;
if ~(xLim(1)<0)
    xLim(1) = -delta;
end
if ~(xLim(2)>0)
    xLim(2) = +delta;
end
xlim(xLim)

yLim = ylim;
if ~(yLim(1)<0)
    yLim(1) = -delta;
end
if ~(yLim(2)>0)
    yLim(2) = +delta;
end
ylim(yLim)

uistack(plot([0 0],ylim,'-k'),'bottom');
uistack(plot(xlim,[0 0],'-k'),'bottom');
grid on
title('after cartNorm')
xlabel('real')
ylabel('imag')
drawnow
suptitle('voxels (repetitions averaged)')

function f = plotPolNormExampleRep(x,y,SVMspace,b)
f = figure('WindowStyle','docked');

%% Polar Normalization
% Plot before
subplot(2,2,1); clear hPP
% polarplot(angle(mean(x,1)),abs(mean(x,1)),'.k'); hold on
polarplot(angle(x(y==1,b)),abs(x(y==1,b)),'.'); hold on
polarplot(angle(x(y==2,b)),abs(x(y==2,b)),'.'); hold on
drawnow
% Normalize
x = polarSpaceNormalization(x,SVMspace);
% Plot after
subplot(2,2,2);
% polarplot(angle(mean(x,1)),abs(mean(x,1)),'.k'); hold on
polarplot(angle(x(y==1,b)),abs(x(y==1,b)),'.'); hold on
polarplot(angle(x(y==2,b)),abs(x(y==2,b)),'.'); hold on
drawnow

%% Cartesian Normalization
% Plot before
subplot(2,2,3);
% hScat = scatter(real(mean(x,1)),imag(mean(x,1)),'ko','filled'); hold on
hScat(1) = scatter(real(x(y==1,b)),imag(x(y==1,b)),'o','filled'); hold on
hScat(2) = scatter(real(x(y==2,b)),imag(x(y==2,b)),'o','filled'); hold on
set(hScat,'MarkerEdgeColor','none');
alpha(hScat,0.5)
ax = gca;
ax.DataAspectRatio = [1 1 1];
ax.PlotBoxAspectRatio = [1 1 1];
xLim = xlim;
delta = abs(diff(xLim)).*0.1;
if ~(xLim(1)<0)
    xLim(1) = -delta;
end
if ~(xLim(2)>0)
    xLim(2) = +delta;
end
xlim(xLim)

yLim = ylim;
if ~(yLim(1)<0)
    yLim(1) = -delta;
end
if ~(yLim(2)>0)
    yLim(2) = +delta;
end
ylim(yLim)

uistack(plot([0 0],ylim,'-k'),'bottom');
uistack(plot(xlim,[0 0],'-k'),'bottom');
grid on
title('before cartNorm')
xlabel('real')
ylabel('imag')
drawnow

% Normalize
x = cartSpaceNormalization(x,SVMspace);

%Plot after
subplot(2,2,4);
hScat(1) = scatter(real(x(y==1,b)),imag(x(y==1,b)),'o','filled'); hold on
hScat(2) = scatter(real(x(y==2,b)),imag(x(y==2,b)),'o','filled'); hold on
set(hScat,'MarkerEdgeColor','none');
alpha(hScat,0.5)
ax = gca;
ax.DataAspectRatio = [1 1 1];
ax.PlotBoxAspectRatio = [1 1 1];
xLim = xlim;
delta = abs(diff(xLim)).*0.1;
if ~(xLim(1)<0)
    xLim(1) = -delta;
end
if ~(xLim(2)>0)
    xLim(2) = +delta;
end
xlim(xLim)

yLim = ylim;
if ~(yLim(1)<0)
    yLim(1) = -delta;
end
if ~(yLim(2)>0)
    yLim(2) = +delta;
end
ylim(yLim)

uistack(plot([0 0],ylim,'-k'),'bottom');
uistack(plot(xlim,[0 0],'-k'),'bottom');
grid on
title('after cartNorm')
xlabel('real')
ylabel('imag')
drawnow
suptitle('voxels (repetitions averaged)')

% function f = plotPolNormExample(x,SVMspace)
% x = permute(x,[2 3 1 4]);
% x = x(:,:,:);
% x = permute(x,[3 4 2 1]);
% f = figure('WindowStyle','docked');
%
% %% Polar Normalization
% subplot(2,2,1); clear hPP
% % polarplot(angle(x(:)),abs(x(:)),'.'); hold on
% for condInd = 1:size(x,3)
%     hPP(condInd) = polarplot(angle(x(:,:,condInd)),abs(x(:,:,condInd)),'o'); hold on
%     hPP(condInd).MarkerFaceColor = hPP(condInd).Color;
%     hPP(condInd).MarkerEdgeColor = 'w';
%     hPP(condInd).MarkerSize = 3.5;
% end
% drawnow
% hPP1 = hPP; clear hPP
% ax1 = gca;
%
% % Normalize
% x = polarSpaceNormalization(x,SVMspace);
%
% % Plot after
% subplot(2,2,2);
% % polarplot(angle(xAfter(:)),abs(xAfter(:)),'.'); hold on
% for condInd = 1:size(x,3)
%     hPP(condInd) = polarplot(angle(x(:,:,condInd)),abs(x(:,:,condInd)),'o'); hold on
%     hPP(condInd).MarkerFaceColor = hPP(condInd).Color;
%     hPP(condInd).MarkerEdgeColor = 'w';
%     hPP(condInd).MarkerSize = 3.5;
% end
% drawnow
% hPP2 = hPP; clear hPP
% ax2 = gca;
%
% ax = ax1;
% ax.ThetaTickLabel = 12-ax.ThetaTick(1:end)/360*12;
% ax.ThetaTickLabel(1,:) = '0 ';
% ax.ThetaAxis.Label.String = {'delay' '(sec)'};
% ax.ThetaAxis.Label.Rotation = 0;
% ax.ThetaAxis.Label.HorizontalAlignment = 'left';
% ax.RAxis.Label.String = 'amp (%BOLD)';
% ax.RAxis.Label.Rotation = 80;
% ax.Title.String = 'before polNorm';
%
% ax = ax2;
% ax.ThetaTickLabel = (-wrapTo180(ax.ThetaTick(1:end))/360*12);
% ax.ThetaTickLabel(1,:) = '0 ';
% % ax.ThetaAxis.Label.String = {'delay' '(sec)'};
% % ax.ThetaAxis.Label.Rotation = 0;
% % ax.ThetaAxis.Label.HorizontalAlignment = 'left';
% ax.Title.String = 'after polNorm';
%
% %% Cartesian Normalization
% subplot(2,2,3); clear hPP
% for condInd = 1:size(x,3)
%     hScat(condInd) = scatter(real(x(:,:,condInd)),imag(x(:,:,condInd)),'o'); hold on
%     hScat(condInd).MarkerFaceColor = hScat(condInd).CData;
%     hScat(condInd).MarkerEdgeColor = 'w';
% %     hScat(condInd).SizeData = 35;
% end
% ax = gca;
% ax.DataAspectRatio = [1 1 1];
% ax.PlotBoxAspectRatio = [1 1 1];
% xLim = xlim;
% delta = abs(diff(xLim)).*0.1;
% if ~(xLim(1)<0)
%     xLim(1) = -delta;
% end
% if ~(xLim(2)>0)
%     xLim(2) = +delta;
% end
% xlim(xLim)
%
% yLim = ylim;
% if ~(yLim(1)<0)
%     yLim(1) = -delta;
% end
% if ~(yLim(2)>0)
%     yLim(2) = +delta;
% end
% ylim(yLim)
%
% uistack(plot([0 0],ylim,'-k'),'bottom');
% uistack(plot(xlim,[0 0],'-k'),'bottom');
% grid on
% title('before cartNorm')
% xlabel('real')
% ylabel('imag')
%
%
% %normalize
% x = cartSpaceNormalization(x,SVMspace);
%
% %after
% subplot(2,2,4); clear hPP
% for condInd = 1:size(x,3)
%     hScat(condInd) = scatter(real(x(:,:,condInd)),imag(x(:,:,condInd)),'o'); hold on
%     hScat(condInd).MarkerFaceColor = hScat(condInd).CData;
%     hScat(condInd).MarkerEdgeColor = 'w';
% %     hScat(condInd).SizeData = 35;
% end
%
% ax = gca;
% ax.DataAspectRatio = [1 1 1];
% ax.PlotBoxAspectRatio = [1 1 1];
% xLim = xlim;
% delta = abs(diff(xLim)).*0.1;
% if ~(xLim(1)<0)
%     xLim(1) = -delta;
% end
% if ~(xLim(2)>0)
%     xLim(2) = +delta;
% end
% xlim(xLim)
%
% yLim = ylim;
% if ~(yLim(1)<0)
%     yLim(1) = -delta;
% end
% if ~(yLim(2)>0)
%     yLim(2) = +delta;
% end
% ylim(yLim)
%
% uistack(plot([0 0],ylim,'-k'),'bottom');
% uistack(plot(xlim,[0 0],'-k'),'bottom');
% grid on
% title('after cartNorm')
% xlabel('real')
% ylabel('imag')
% drawnow
%
% suptitle('repetitions (most active voxel)')



function printRes(resSess,resSubj,resGroup)
disp(['Session results:'])
disp([' -acc (one-sided p)'])
disp(resSess.acc)
disp(resSess.acc_p)
disp([' -distT (one-sided p)'])
disp(resSess.distT)
disp(resSess.distT_p)

disp(['Group results:'])
disp(['hit    =' num2str(resGroup.hit) '/' num2str(resGroup.nObs)])
disp([' -fixedEffect'])
disp(['  acc    =' num2str(resGroup.acc*100,'%0.2f%%') '; p=' num2str(resGroup.acc_p,'%0.3f') '; thresh=' num2str(resGroup.acc_thresh,'%0.2f')])
disp(['  auc    =' num2str(resGroup.auc*100,'%0.2f') '; 90%CI=' num2str([resGroup.auc_CI5 resGroup.auc_CI95],'%0.3f ')])
disp(['  distT  =' num2str(resGroup.distT,'%0.2f') '; p=' num2str(resGroup.distT_p,'%0.3f ')])
disp([' -randomEffect'])
disp(['  acc=' num2str(mean(resSubj.acc)*100,'%0.2f%%')])
disp(['  -student'])
disp(['   T=' num2str(resGroup.acc_T,'%0.2f') '; P=' num2str(resGroup.acc_P,'%0.3f')])
disp(['  -wilcoxon'])
disp(['   sRank=' num2str(resGroup.acc_wilcoxonSignedrank,'%0.2f') '; P=' num2str(resGroup.acc_wilcoxonP,'%0.3f')])

function printRes2(resGroup)
disp(['FixedEffect'])
if isfield(resGroup,'auc')
    disp([' -auc    =' num2str(resGroup.auc*100,'%0.2f') '; 90%CI=' num2str([resGroup.auc_CI5 resGroup.auc_CI95],'%0.3f ')])
end
if isfield(resGroup,'acc')
    disp([' -acc    =' num2str(resGroup.acc*100,'%0.2f') '; 90%CI=' num2str([resGroup.acc_CI5 resGroup.acc_CI95],'%0.3f ')])
end
disp(['RandomEffect'])
if isfield(resGroup,'distT_T')
    disp([' -student on distT'])
    disp(['  T=' num2str(resGroup.distT_T,'%0.2f') '; ones-sided P=' num2str(resGroup.distT_P,'%0.3f')])
end
if isfield(resGroup,'auc_T')
    disp([' -student on auc'])
    disp(['  T=' num2str(resGroup.auc_T,'%0.2f') '; ones-sided P=' num2str(resGroup.auc_P,'%0.3f')])
end
if isfield(resGroup,'acc_T')
    disp([' -student on acc'])
    disp(['  T=' num2str(resGroup.acc_T,'%0.2f') '; ones-sided P=' num2str(resGroup.acc_P,'%0.3f')])
end
if isfield(resGroup,'auc_wilcoxonSignedrank')
    disp([' -wilcoxon on auc'])
    disp(['  sRank=' num2str(resGroup.auc_wilcoxonSignedrank,'%0.2f') '; ones-sided P=' num2str(resGroup.auc_wilcoxonP,'%0.3f')])
end
if isfield(resGroup,'acc_wilcoxonSignedrank')
    disp([' -wilcoxon on acc'])
    disp(['  sRank=' num2str(resGroup.acc_wilcoxonSignedrank,'%0.2f') '; ones-sided P=' num2str(resGroup.acc_wilcoxonP,'%0.3f')])
end



 1
 2
 3
 4
 5
 6
 7
 8
 9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
78
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97

	

function res = perfMetric(y,yHat,k,liteFlag)
warning('off','stats:perfcurve:SubSampleWithMissingClasses')
averageWR = 1;
if ~exist('liteFlag','var')
    liteFlag = false;
end
if ~exist('y','var')
    res = struct(...
        'y',[],...
        'yHat',[],...
        'nObs',[],...
        'hit',[],...
        'acc',[],...
        'acc_CI5',[],...
        'acc_CI95',[],...
        'acc_thresh',[],...
        'acc_p',[],...);
        'auc',[],...
        'auc_CI5',[],...
        'auc_CI95',[],...
        'distT',[],...
        'distT_p',[],...
        'distMean',[],...
        'distStd',[]);
    return
end

if averageWR
    nRun = length(unique(k))*length(unique(y));
    y = mean(reshape(y,[length(y)/nRun nRun]),1)';
    if nRun~=size(yHat,1)
        error('double-check that')
    end
    yHat = mean(reshape(yHat,[size(yHat,1)/nRun nRun size(yHat,3) size(yHat,4)]),1);
    yHat = permute(yHat,[2 3 4 1]); % run x t x wSect
end

res.y = {y};
res.nObs = length(y);
res.yHat = {yHat};

% Contitnue only if we have a pair of conditions
if length(unique(y))>2
    resTmp = perfMetric;
    fieldList = fields(res);
    for fieldInd = 1:length(fieldList)
        resTmp.(fieldList{fieldInd}) = res.(fieldList{fieldInd});
    end
    res = resTmp;
    return
end

% Lite version
if liteFlag
    res = perfMetric;
    for tInd = 1:size(yHat,2)
        % auc
        [~,~,~,auc] = perfcurve(y,yHat(:,tInd),1);
        res.auc(:,tInd) = auc(1);
    end
    return
end


% acc
if size(yHat,2)==1
    res.hit = sum((yHat<0)+1==y);
    res.acc = res.hit./res.nObs;
    [~,pci] = binofit(res.hit,res.nObs,0.1);
    res.acc_CI5 = pci(1);
    res.acc_CI95 = pci(2);
    [~,pci] = binofit(res.nObs/2,res.nObs,0.1);
    res.acc_thresh = pci(2);
    res.acc_p = binocdf(res.hit,res.nObs,0.5,'upper');
else
    res.hit = [];
    res.acc = [];
    res.acc_CI5 = [];
    res.acc_CI95 = [];
    res.acc_thresh = [];
    res.acc_p = [];
end
for tInd = 1:size(yHat,2)
    % auc
    [~,~,~,auc] = perfcurve(y,yHat(:,tInd),1,'NBOOT',2^10);
    res.auc(:,tInd) = auc(1);
    res.auc_CI5(:,tInd) = auc(2);
    res.auc_CI95(:,tInd) = auc(3);
    % distT
    [~,P,~,STATS] = ttest(yHat(y==1,tInd),yHat(y==2,tInd));
    res.distT(:,tInd) = STATS.tstat;
    res.distT_p(:,tInd) = P;
end
% dist
dist = yHat(y==1,:) - yHat(y==2,:);
res.distMean = mean(dist,1);
res.distStd = std(dist,[],1);


function [resSess,resSubj,resGroup] = summarizePerf(res_sess)
allField = fields(res_sess);
for i = 1:length(allField)
    if isnumeric(res_sess(1).(allField{i}))
        tLength = length(res_sess(1).(allField{i}));
        if tLength==1
            resSess.(allField{i}) = nan(size(res_sess));
        elseif tLength==0
            resSess.(allField{i}) = [];
        else % for the special case where we compute stats at each time delay relative to stimulus onset
            resSess.(allField{i}) = nan([tLength size(res_sess)]);
        end
        if ~isempty(res_sess(1).(allField{i}))
            resSess.(allField{i})(:) = [res_sess.(allField{i})];
        end
        if length(size(resSess.(allField{i})))>2
            resSess.(allField{i}) = permute(resSess.(allField{i}),[2 3 1]);
        end
    elseif iscell(res_sess(1).(allField{i}))
        resSess.(allField{i}) = cell(size(res_sess));
        resSess.(allField{i})(:) = [res_sess.(allField{i})];
    elseif ischar(res_sess(1).(allField{i}))
        resSess.(allField{i}) = cell(size(res_sess));
        resSess.(allField{i})(:) = {res_sess.(allField{i})};
    elseif isstruct(res_sess(1).(allField{i}))
        resSess.(allField{i}) = cell(size(res_sess));
        resSess.(allField{i})(:) = {res_sess.(allField{i})};
    else
        error('code that')
    end
end
resSess.info = 'subj x sess';
resSess.distT_fdr = resSess.distT_p;
resSess.distT_fdr(:) = mafdr(resSess.distT_p(:),'BHFDR',true);
% resSess = orderfields(resSess,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 23 16 17 18 19 20 21 22]);

% Continue only if having a pair of condtions
if isempty(resSess.acc)
    resSubj  = [];
    resGroup = [];
    return
end

resSubj.y = cell(size(resSess.y,1),1);
resSubj.yHat = cell(size(resSess.yHat,1),1);
for subjInd = 1:size(resSess.y,1)
    resSubj.y{subjInd} = cell2mat(resSess.y(subjInd,:)');
    resSubj.yHat{subjInd} = cell2mat(resSess.yHat(subjInd,:)');
end
resSubj.nObs = sum(resSess.nObs,2);
if ~isempty(resSess.hit)
    resSubj.hit = sum(resSess.hit,2);
    resSubj.acc = resSubj.hit./resSubj.nObs;
    [~,pci] = binofit(resSubj.hit,resSubj.nObs,0.1);
    resSubj.acc_CI5 = pci(:,1);
    resSubj.acc_CI95 = pci(:,2);
    [~,pci] = binofit(resSubj.nObs/2,resSubj.nObs,0.1);
    resSubj.acc_thresh = pci(:,2);
    resSubj.acc_p = binocdf(resSubj.hit,resSubj.nObs,0.5,'upper');
    resSubj.acc_fdr = mafdr(resSubj.acc_p,'BHFDR',true);
else
    resSubj.hit = [];
    resSubj.acc = [];
    resSubj.acc_CI5 = [];
    resSubj.acc_CI95 = [];
    resSubj.acc_thresh = [];
    resSubj.acc_p = [];
    resSubj.acc_fdr = [];
end
resSubj.auc = nan([size(resSubj.y) size(resSubj.yHat{1},2)]);
resSubj.distT = nan([size(resSubj.y) size(resSubj.yHat{1},2)]);
resSubj.distT_p = nan([size(resSubj.y) size(resSubj.yHat{1},2)]);
resSubj.distMean = nan([size(resSubj.y) size(resSubj.yHat{1},2)]);
resSubj.distStd = nan([size(resSubj.y) size(resSubj.yHat{1},2)]);
for subjInd = 1:size(resSubj.y,1)
    for tInd = 1:size(resSubj.yHat{subjInd},2)
        [~,~,~,auc] = perfcurve(resSubj.y{subjInd},resSubj.yHat{subjInd}(:,tInd),1);
        resSubj.auc(subjInd,:,tInd) = auc(1);
        resSubj.auc_CI5(subjInd,:,tInd) = nan;
        resSubj.auc_CI95(subjInd,:,tInd) = nan;
        %     [~,~,~,auc] = perfcurve(resSubj.y{subjInd},resSubj.yHat{subjInd},1,'NBOOT',2^10);
        %     resSubj.auc(subjInd) = auc(1);
        %     resSubj.auc_CI5(subjInd) = auc(2);
        %     resSubj.auc_CI95(subjInd) = auc(3);
        
        [~,P,~,STATS] = ttest(resSubj.yHat{subjInd}(resSubj.y{subjInd}==1,tInd),resSubj.yHat{subjInd}(resSubj.y{subjInd}==2,tInd));
        resSubj.distT(subjInd,:,tInd) = STATS.tstat;
        resSubj.distT_p(subjInd,:,tInd) = P;
        
%         resSubj.distMean(subjInd,:,tInd) = mean(resSess.distMean(subjInd,:,tInd),2);
%         resSubj.distStd(subjInd,:,tInd) = std(resSess.distMean(subjInd,:,tInd),[],2);
    end
    resSubj.distMean(subjInd,:,:) = mean(resSess.distMean(subjInd,:,:),2);
    resSubj.distStd(subjInd,:,:) = std(resSess.distMean(subjInd,:,:),[],2);
end
resSubj.nVoxOrig = round(median(resSess.nVoxOrig,2));
resSubj.nVox = round(median(resSess.nVox,2));
resSubj.chanSpace = resSess.chanSpace(:,1);
resSubj.complexSpace = resSess.complexSpace(:,1);
resSubj.svmKernel = resSess.svmKernel(:,1);

resGroup.y = cell2mat(resSubj.y);
resGroup.yHat = cell2mat(resSubj.yHat);
resGroup.nObs = sum(resSubj.nObs,1);
if ~isempty(resSubj.hit)
    resGroup.hit = sum(resSubj.hit,1);
    
    resGroup.acc = resGroup.hit./resGroup.nObs;
    [~,pci] = binofit(resGroup.hit,resGroup.nObs,0.1);
    resGroup.acc_CI5 = pci(:,1);
    resGroup.acc_CI95 = pci(:,2);
    [~,pci] = binofit(resGroup.nObs/2,resGroup.nObs,0.1);
    resGroup.acc_thresh = pci(:,2);
    resGroup.acc_p = binocdf(resGroup.hit,resGroup.nObs,0.5,'upper');
    [~,P,~,STATS] = ttest(resSubj.acc,0.5,'tail','right');
    resGroup.acc_T = STATS.tstat;
    resGroup.acc_P = P;
    [P,~,STATS] = signrank(resSubj.acc,0.5,'tail','right');
    resGroup.acc_wilcoxonSignedrank = STATS.signedrank;
    resGroup.acc_wilcoxonP = P;
else
    resGroup.hit = [];
    
    resGroup.acc = [];
    resGroup.acc_CI5 = [];
    resGroup.acc_CI95 = [];
    resGroup.acc_thresh = [];
    resGroup.acc_p = [];
    resGroup.acc_T = [];
    resGroup.acc_P = [];
    resGroup.acc_wilcoxonSignedrank = [];
    resGroup.acc_wilcoxonP = [];
end

resGroup.auc = nan(1,1,size(resGroup.yHat,2));
resGroup.auc_CI5 = nan(1,1,size(resGroup.yHat,2));
resGroup.auc_CI95 = nan(1,1,size(resGroup.yHat,2));
resGroup.auc_T = nan(1,1,size(resGroup.yHat,2));
resGroup.auc_P = nan(1,1,size(resGroup.yHat,2));
resGroup.auc_wilcoxonSignedrank = nan(1,1,size(resGroup.yHat,2));
resGroup.auc_wilcoxonP = nan(1,1,size(resGroup.yHat,2));
resGroup.distT = nan(1,1,size(resGroup.yHat,2));
resGroup.distT_p = nan(1,1,size(resGroup.yHat,2));
resGroup.distT_T = nan(1,1,size(resGroup.yHat,2));
resGroup.distT_P = nan(1,1,size(resGroup.yHat,2));
resGroup.distMean_T = nan(1,1,size(resGroup.yHat,2));
resGroup.distMean_P = nan(1,1,size(resGroup.yHat,2));
for tInd = 1:size(resGroup.yHat,2)
    [~,~,~,auc] = perfcurve(resGroup.y,resGroup.yHat(:,tInd),1,'NBOOT',2^10);
    resGroup.auc(tInd) = auc(1);
    resGroup.auc_CI5(tInd) = auc(2);
    resGroup.auc_CI95(tInd) = auc(3);
    [~,P,~,STATS] = ttest(resSubj.auc(:,:,tInd),0.5,'tail','right');
    resGroup.auc_T(tInd) = STATS.tstat;
    resGroup.auc_P(tInd) = P;
    [P,~,STATS] = signrank(resSubj.auc(:,:,tInd),0.5,'tail','right');
    resGroup.auc_wilcoxonSignedrank(tInd) = STATS.signedrank;
    resGroup.auc_wilcoxonP(tInd) = P;

    [~,P,~,STATS] = ttest(resGroup.yHat(resGroup.y==1,tInd),resGroup.yHat(resGroup.y==2,tInd),'tail','right');
    resGroup.distT(tInd) = STATS.tstat;
    resGroup.distT_p(tInd) = P;
    [~,P,~,STATS] = ttest(resSubj.distT(:,:,tInd),0,'tail','right');
    resGroup.distT_T(tInd) = STATS.tstat;
    resGroup.distT_P(tInd) = P;
    [~,P,~,STATS] = ttest(resSubj.distMean(:,:,tInd),0,'tail','right');
    resGroup.distMean_T(tInd) = STATS.tstat;
    resGroup.distMean_P(tInd) = P;
end
resGroup.nVoxOrig = round(median(resSubj.nVoxOrig,1));
resGroup.nVox = round(median(resSubj.nVox,1));
resGroup.chanSpace = resSubj.chanSpace{1};

