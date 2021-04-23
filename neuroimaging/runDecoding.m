function [resBS,resWS] = runDecoding(p,verbose,nPerm,figOption)
if ~exist('verbose','var')
    verbose = 1;
end
if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end
if ~isfield(p,'svmSpace') || isempty(p.svmSpace)
    p.svmSpace = 'cart'; % 'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'
    % 'hr' 'hrNoAmp' 'cart' 'cartNoAmp' cartNoAmp_HT 'cartReal', 'cartImag', 'pol', 'polMag' 'polMag_T' or 'polDelay'
end
if isstruct(p.svmSpace)
    doPerm = 1;
    res = p.svmSpace;
    p.svmSpace = res.summary.svmSpace;
    if ~exist('nPerm','var') || isempty(nPerm)
        nPerm = 100;
    end
else
    doPerm = 0;
end
if isfield(p,'dataType') && ~isempty(p.dataType)
    dataType = p.dataType;
else
    dataType = 'sin';
end
 


%% Define paths
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'};
if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
        funPath = fullfile(repoPath,'C-derived\DecodingHR\fun');
            inDir  = 'd';
%make sure everything is forward slash for mac, linux pc compatibility
for tmp = {'repoPath' 'funPath' 'inDir'}
    eval([char(tmp) '(strfind(' char(tmp) ',''\''))=''/'';']);
end
clear tmp


if doPerm
    error('double-check that')
    filename = fullfile(pwd,mfilename);
    filename = fullfile(filename,[p.svmSpace '_' num2str(nPerm) 'perm']);
    if exist([filename '.mat'],'file')
        load(filename)
        if verbose
            disp('permutation found on disk, skipping')
            disp('Group results:')
            disp(['  hit    =' num2str(res.summary.hit) '/' num2str(res.summary.nObs)])
            disp(['  acc    =' num2str(res.summary.acc*100,'%0.2f%%')])
            disp(' permutation test stats')
            disp(['  thresh =' num2str(res.perm.summary.accThresh*100,'%0.2f%%')])
            disp(['  p      =' num2str(res.perm.summary.p,'%0.3f')])
        end
        return
    else
        disp(['running ' num2str(nPerm) ' permutations']);
    end
end


%% Load data
dAll = cell(size(subjList,1),1);
for subjInd = 1:size(subjList,2)
    curFile = fullfile(funPath,inDir,[subjList{subjInd} '.mat']);
    if verbose; disp(['loading: ' curFile]); end
    load(curFile,'res');
    dAll{subjInd} = res;
end
d = dAll; clear dAll
sessList = fields(d{1});
% Load feature slection
load(fullfile(funPath,inDir,'featSel.mat'),'featSel');
if verbose
    disp('---');
    disp(['SVM space: ' p.svmSpace '; on ' dataType]);
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

%% Define feature selection
featSelSteps_labelList = featSel{1,1}.featSeq.featSelList;
for i = 1:length(featSelSteps_labelList)
    tmp = strsplit(featSelSteps_labelList{i},':');
    featSelSteps_labelList(i) = tmp(1);
end
featSelConds_labelList = featSel{1,1}.featSeq.condPairList;
for i = 1:length(featSelConds_labelList)
    featSelConds_labelList{i} = num2str(featSelConds_labelList{i});
end
switch p.featSel.global.method
    case 'allData'
        % Non-condition-specific featSel steps
        ind_nSpecFeatSel = ismember(featSelSteps_labelList,{'retinoFov' 'act' 'respVecSig' 'vein' 'respVecDiff'});
        % which condition set
        ind_nSpecFeatSelCond = [1 2 3];
        % Condition-specific fetSel steps
        ind_specFeatSel = ismember(featSelSteps_labelList,{''});
    case 'custom1'
        % Non-condition-specific featSel steps
        ind_nSpecFeatSel = ismember(featSelSteps_labelList,{'retinoFov' 'act' 'respVecSig' 'vein'});
        % which condition set
        ind_nSpecFeatSelCond = [1 2 3];
        % Condition-specific fetSel steps
        ind_specFeatSel = ismember(featSelSteps_labelList,{'respVecDiff'});
        % which condition set
        switch p.condPair
            case 'grat1VSgrat2'
                ind_specFeatSelCond = [1 2];
            case 'grat1VSplaid'
                ind_specFeatSelCond = [1 3];
            case 'grat2VSplaid'
                ind_specFeatSelCond = [2 3];
            otherwise
                error('X')
        end
    case 'custom2'
        % Non-condition-specific featSel steps
        ind_nSpecFeatSel = ismember(featSelSteps_labelList,{'retinoFov' 'vein'});
        % which condition set
        ind_nSpecFeatSelCond = [1 2 3];
        % Condition-specific fetSel steps
        ind_specFeatSel = ismember(featSelSteps_labelList,{'act' 'respVecSig' 'respVecDiff'});
        % which condition set
        switch p.condPair
            case 'grat1VSgrat2'
                ind_specFeatSelCond = [1 2];
            case 'grat1VSplaid'
                ind_specFeatSelCond = [1 3];
            case 'grat2VSplaid'
                ind_specFeatSelCond = [2 3];
            otherwise
                error('X')
        end
    otherwise
        error('X')
end
ind_nSpecFeatSelCond = squeeze(ismember(featSelConds_labelList,num2str(ind_nSpecFeatSelCond)));
ind_specFeatSelCond = squeeze(ismember(featSelConds_labelList,num2str(ind_specFeatSelCond)));

for subjInd = 1:size(d,1)
    for sessInd = 1:size(d,2)
        featSel{subjInd,sessInd}.indIn = ...
            all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_nSpecFeatSel,ind_nSpecFeatSelCond),2)...
            & all(featSel{subjInd,sessInd}.featSeq.featIndIn(:,ind_specFeatSel,ind_specFeatSelCond),2);
    end
end


%% Example plot of trigonometric (polar) representation
subjInd = p.figOption.subjInd;
sessInd = p.figOption.sessInd;
if p.figOption.verbose>1
    f = plotNorm(d{subjInd,sessInd},p,featSel{subjInd,sessInd});
    if figOption.save
        error('code that')
        filename = fullfile(pwd,mfilename);
        if ~exist(filename,'dir'); mkdir(filename); end
        filename = fullfile(filename,p.svmSpace);
        f.Color = 'none';
        set(findobj(f.Children,'type','Axes'),'color','none')
        set(findobj(f.Children,'type','PolarAxes'),'color','none')
        saveas(f,[filename '.svg']); if verbose; disp([filename '.svg']); end
        f.Color = 'w';
        set(findobj(f.Children,'type','Axes'),'color','w')
        set(findobj(f.Children,'type','PolarAxes'),'color','w')
        saveas(f,[filename '.fig']); if verbose; disp([filename '.fig']); end
        saveas(f,[filename '.jpg']); if verbose; disp([filename '.jpg']); end
    end
end

%% Within-session SVM cross-validation (with cross-session feature selection)
svmModelK = cell(size(d));
nrmlzK = cell(size(d));
yHatK = cell(size(d));
yHatK_tr = cell(size(d));
Y = cell(size(d));
K = cell(size(d));
nOrig = nan(size(d));
n = nan(size(d));

disp('Within-session SVM cross-validation')
disp('with between-session feature selection')
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
    
    if doPerm
        error('code that')
        disp(['for sess ' num2str(i) ' of ' num2str(numel(d))])
        tic
    end
    
    switch dataType
        case 'sin'
            [X,y,k] = getXYK(d{subjInd,sessInd},p);
        case 'waveTrialSparse' % do not average and use 1 tPts out of 12 in each stimulus cycle
            error('code that')
            [X,y,k,~] = getXYK_wave(d{subjInd,sessInd},SVMspace,'trialSparse');
        case 'waveTrialSparseCat2' % do not average and use 1 tPts out of 12 in each stimulus cycle, concatenating each tPts in the feature (vox) dimension
            error('code that')
            [X,y,k,~] = getXYK_wave(d{subjInd,sessInd},SVMspace,'trialSparseCat2');
        case 'waveTrialSparseRep' % do not average and use 1 tPts out of 12 in each stimulus cycle, repeating svm training for each time points and averaging only the end model
            error('code that')
            [X,y,k,t] = getXYK_wave(d{subjInd,sessInd},SVMspace,'full');
            tmp = repmat(1:12,[length(t)/12 1])';
            t = tmp(:);
        case 'waveRun' % average within runs
            error('code that')
            [X,y,k,~] = getXYK_wave(d{subjInd,sessInd},SVMspace,'run');
        case 'wave' % average within trials
            error('code that')
            [X,y,k,~] = getXYK_wave(d{subjInd,sessInd},SVMspace,'trial');
        case 'waveFull' % do not average and use all tPts
            error('code that')
            [X,y,k,~] = getXYK_wave(d{subjInd,sessInd},SVMspace,'full');
        otherwise
            error('X')
    end
    
    % cross-session feature selection
    featSelInd = featSel{subjInd,sessIndCross}.indIn;
    if strcmp(dataType,'waveTrialSparseCat2')
        featSelInd = repmat(featSelInd,[1 12]);
    end
    x = X(:,featSelInd);
    
    % k-fold training
    if ~strcmp(dataType,'waveTrialSparseRep')
        [svmModelK{i},nrmlzK{i}] = SVMtrain(y,x,p,k);
    else
        [svmModelK{i},nrmlzK{i}] = SVMtrain(y,x,p,k,t);
    end
    % k-fold testing
    [yHatK{i},yHatK_tr{i}] = SVMtest(y,x,svmModelK{i},nrmlzK{i},k);
    yHatK{i} = nanmean(yHatK{i},2);
    % output
    Y{i} = y;
    K{i} = k;
    nOrig(i) = length(featSelInd);
    n(i) = nnz(featSelInd);
end


%% Within-session SVM training (and feature selection) + Cross-session SVM testing
disp('Between-session SVM cross-validation')
disp('with feature selection from the training session')
svmModel = cell(size(d));
nrmlz = cell(size(d));
yHat_tr = cell(size(d));
yHat = cell(size(d));
% Training
for i = 1:numel(d)
    if verbose
        disp(['-training sess ' num2str(i) '/' num2str(numel(d))])
    end
    [subjInd,trainInd] = ind2sub(size(d),i);
    if doPerm
        error('code that')
        disp(['for sess ' num2str(i) ' of ' num2str(numel(d))])
        tic
    end
    % get training data
    switch dataType
        case 'sin'
            [X,y,~] = getXYK(d{subjInd,trainInd},p);
        case 'waveTrialSparse' % do not average and use 1 tPts out of 12 in each stimulus cycle
            error('code that')
            [X,y,~,~] = getXYK_wave(d{subjInd,trainInd},SVMspace,'trialSparse');
        case 'waveTrialSparseCat2' % do not average and use 1 tPts out of 12 in each stimulus cycle, concatenating each tPts in the feature (vox) dimension
            error('code that')
            [X,y,~,~] = getXYK_wave(d{subjInd,trainInd},SVMspace,'trialSparseCat2');
        case 'waveTrialSparseRep' % do not average and use 1 tPts out of 12 in each stimulus cycle, repeating svm training for each time points and averaging only the end model
            error('code that')
            [X,y,~,t] = getXYK_wave(d{subjInd,trainInd},SVMspace,'full');
            tmp = repmat(1:12,[length(t)/12 1])';
            t = tmp(:);
        case 'waveRun' % average within runs
            error('code that')
            [X,y,~,~] = getXYK_wave(d{subjInd,trainInd},SVMspace,'run');
        case 'wave' % average within trials
            error('code that')
            [X,y,~,~] = getXYK_wave(d{subjInd,trainInd},SVMspace,'trial');
        case 'waveFull' % do not average and use all tPts
            error('code that')
            [X,y,~,~] = getXYK_wave(d{subjInd,trainInd},SVMspace,'full');
        otherwise
            error('X')
    end
    % same-session feature-selection
    featSelInd = featSel{subjInd,trainInd}.indIn;
    if strcmp(dataType,'waveTrialSparseCat2')
        featSelInd = repmat(featSelInd,[1 12]);
    end
    x = X(:,featSelInd);
    % training
    if ~strcmp(dataType,'waveTrialSparseRep')
        [svmModel{subjInd,trainInd},nrmlz{subjInd,trainInd}] = SVMtrain(y,x,p);
    else
        [svmModel{subjInd,trainInd},nrmlz{subjInd,trainInd}] = SVMtrain(y,x,p,[],t);
    end
    % same session testing (not crossvalidated)
    [yHat_tr{subjInd,trainInd},~] = SVMtest(y,x,svmModel{subjInd,trainInd},nrmlz{subjInd,trainInd});
end
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
    switch dataType
        case 'sin'
            [X,y,~] = getXYK(d{subjInd,testInd},p);
        case 'waveTrialSparse' % do not average and use 1 tPts out of 12 in each stimulus cycle
            error('code that')
            [X,y,~,~] = getXYK_wave(d{subjInd,testInd},SVMspace,'trialSparse');
        case 'waveTrialSparseCat2' % do not average and use 1 tPts out of 12 in each stimulus cycle
            error('code that')
            [X,y,~,~] = getXYK_wave(d{subjInd,testInd},SVMspace,'trialSparseCat2');
        case 'waveTrialSparseRep' % do not average and use 1 tPts out of 12 in each stimulus cycle, repeating svm training for each time points and averaging only the end model
            error('code that')
            [X,y,~,t] = getXYK_wave(d{subjInd,testInd},SVMspace,'full');
            tmp = repmat(1:12,[length(t)/12 1])';
            t = tmp(:);
        case 'waveRun' % average within runs
            error('code that')
            [X,y,~,~] = getXYK_wave(d{subjInd,testInd},SVMspace,'run');
        case 'wave' % average within trials
            error('code that')
            [X,y,~,~] = getXYK_wave(d{subjInd,testInd},SVMspace,'trial');
        case 'waveFull' % do not average and use all tPts
            error('code that')
            [X,y,~,~] = getXYK_wave(d{subjInd,testInd},SVMspace,'full');
        otherwise
            error('X')
    end
    % cross-session feature selection
    featSelInd = featSel{subjInd,trainInd}.indIn;
    if strcmp(dataType,'waveTrialSparseCat2')
        featSelInd = repmat(featSelInd,[1 12]);
    end
    x = X(:,featSelInd);
    % cross-session SVM testing
    [yHat{subjInd,testInd},~] = SVMtest(y,x,svmModel{subjInd,trainInd});
end

%% Compute performance metrics
% Between-sess
resBS_sess = repmat(perfMetric,[size(d,1) size(d,2)]);
for i = 1:numel(d)
    resBS_sess(i) = perfMetric(Y{i},yHat{i},K{i});
end
acc_fdr = mafdr([resBS_sess.acc_p],'BHFDR',true);
for i = 1:numel(d)
    resBS_sess(i).acc_fdr = acc_fdr(i);
    resBS_sess(i).nVoxOrig = size(d{i}.sin,2);
    resBS_sess(i).nVox = nnz(featSel{i}.indIn);
    resBS_sess(i).svmSpace = p.svmSpace;
    resBS_sess(i).norm = p.norm;
end
resBS_sess = orderfields(resBS_sess,[1 2 3 4 5 6 7 8 9 15 10 11 12 13 14 16 17 18 19]);
                                    
% Within-sess
resWS_sess = repmat(perfMetric,[size(d,1) size(d,2)]);
for i = 1:numel(d)
    resWS_sess(i) = perfMetric(Y{i},yHatK{i},K{i});
end
acc_fdr = mafdr([resWS_sess.acc_p],'BHFDR',true);
for i = 1:numel(d)
    resWS_sess(i).acc_fdr = acc_fdr(i);
    resWS_sess(i).nVoxOrig = size(d{i}.sin,2);
    resWS_sess(i).nVox = nnz(featSel{i}.indIn);
    resWS_sess(i).svmSpace = p.svmSpace;
    resWS_sess(i).norm = p.norm;
end
resWS_sess = orderfields(resWS_sess,[1 2 3 4 5 6 7 8 9 15 10 11 12 13 14 16 17 18 19]);

%% Summarize group performances
[resBSsess,resBSsubj,resBSgroup] = summarizePerf(resBS_sess);
[resWSsess,resWSsubj,resWSgroup] = summarizePerf(resWS_sess);
fieldList = fields(resBSsubj);
for i = 1:length(fieldList)
    if isnumeric(resBSsubj.(fieldList{i}))
        resSubj.(fieldList{i}) = mean(cat(3,resWSsubj.(fieldList{i}),resBSsubj.(fieldList{i})),3);
    end
end
[~,P,~,STATS] = ttest(resSubj.acc,0.5,'tail','right');
resGroup.acc_T = STATS.tstat;
resGroup.acc_P = P;
[~,P,~,STATS] = ttest(resSubj.auc,0.5,'tail','right');
resGroup.auc_T = STATS.tstat;
resGroup.auc_P = P;
[P,~,STATS] = signrank(resSubj.acc,0.5,'tail','right');
resGroup.acc_wilcoxonSignedrank = STATS.signedrank;
resGroup.acc_wilcoxonP = P;
[P,~,STATS] = signrank(resSubj.auc,0.5,'tail','right');
resGroup.auc_wilcoxonSignedrank = STATS.signedrank;
resGroup.auc_wilcoxonP = P;

%% Print info and output
if verbose
    disp('-----------------')
    disp('*Within-session*')
    printRes2(resWSgroup)
    disp(' ')
    disp('*Within+Between*')
    printRes2(resGroup)
    disp(' ')
    disp('*Between-session*')
    printRes2(resBSgroup)
    disp('-----------------')
end
resBS.sess = resBSsess;
resBS.sess.subjList = subjList;
resBS.subj = resBSsubj;
resBS.subj.subjList = subjList;
resBS.group = resBSgroup;
resWS.sess = resWSsess;
resWS.sess.subjList = subjList;
resWS.subj = resWSsubj;
resWS.subj.subjList = subjList;
resWS.group = resWSgroup;

resBS.sess.nVoxOrig = nOrig;
resWS.sess.nVoxOrig = n;
resBS.sess.nVox = n;
resWS.sess.nVox = n;

%% Plot between-session vs within-session
if verbose
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
    
    figure('WindowStyle','docked');
    compareRes(resBS,resWS)
    xlabel('between-session')
    ylabel('within-session')
    textLine = {};
    if p.figOption.verbose>1
        textLine = [textLine {[p.svmSpace '; ' p.condPair '; ' dataType]}];
    else
        textLine = [textLine {[p.svmSpace '; ' p.condPair]}];
    end
    if p.figOption.verbose>1
        textLine = [textLine {strjoin(featSel{figOption.subj,1}.featSeq.featSelList,'; ')}];
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
switch p.svmSpace
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
normTr.svmSpace = '';
normTr.pol = struct('rhoScale',[],'thetaShift',[]);
normTr.svm = struct('scale',[],'shift',[]);
normTr = repmat(normTr,[1 1 length(kList)]);
svmModel = repmat(struct('k',[],'svmSpace',[],'norm',[],'w',[],'b',[],'Paramters',[],'info',[]),[1 1 length(kList)]);
for kInd = 1:length(kList)
    % Split train and test
    te = K==kList(kInd);
    y = Y(~te,:);
    x = X(~te,:);
    
    % Normalize
    [x,normTr(:,:,kInd).pol] = polarSpaceNormalization(x,p.svmSpace);
    [x,normTr(:,:,kInd).svm] = cartSpaceNormalization(x,p.svmSpace,[],p.norm.doCartSpaceScale);
    [x,~,~] = complex2svm(x,p.svmSpace);
    normTr(:,:,kInd).k = kList(kInd);
    normTr(:,:,kInd).svmSpace = p.svmSpace;
    
    % Train
    model = svmtrain(y,x,'-t 0 -q');
    w = model.sv_coef'*model.SVs;
    b = model.rho;
    toUnitNorm = norm(w);
    w = w./toUnitNorm;
    b = b./toUnitNorm;
    
    % Store
    svmModel(:,:,kInd).k = kList(kInd);
    svmModel(:,:,kInd).w = w;
    svmModel(:,:,kInd).b = b;
    svmModel(:,:,kInd).Paramters = model.Parameters;
    svmModel(:,:,kInd).info = '   yHat = x*w''-b   ';
    svmModel(:,:,kInd).svmSpace = p.svmSpace;
    svmModel(:,:,kInd).norm = p.norm;
end

function [yHat,yHatTr] = SVMtest(Y,X,svmModel,nrmlz,K)
if ~exist('K','var') || isempty(K) || all(K(:)==0)
    K = ones(size(Y)); % no crossvalidation
end
kList = [svmModel.k];

yHat = nan([size(Y,1) length(kList)]);
yHatTr = nan([size(Y,1) length(kList)]);
for kInd = 1:length(kList)
    x = X;
    SVMspace = svmModel(kInd).svmSpace;
    norm = svmModel(kInd).norm;
    % Split train and test
    if kList(kInd)==0
        te = true(size(Y));
    else
        te = K==kList(kInd);
    end
    
    % Normalize
    if ~exist('nrmlz','var') || isempty(nrmlz)
        [x,~] = polarSpaceNormalization(x,svmModel(kInd).svmSpace);
        [x,~] = cartSpaceNormalization(x,svmModel(kInd).svmSpace,[],norm.doCartSpaceScale);
    else
        [x,~] = polarSpaceNormalization(x,nrmlz(kInd),te);
        [x,~] = cartSpaceNormalization(x,nrmlz(kInd),te);
    end
    [x,~,~] = complex2svm(x,SVMspace);
    
    w = svmModel(kInd).w;
    b = svmModel(kInd).b;
    yHat(te,kInd) = x(te,:)*w'-b;
    yHatTr(~te,kInd) = x(~te,:)*w'-b;
end
    


% function [yTe,d,yHatTe] = xValSVM(x,y,k,SVMspace)
% X = x;
% kList = unique(k);
% % yTr = nan(length(y),length(kList));
% yTe = nan(length(y),1);
% % yHatTr = nan(length(y),length(kList));
% yHatTe = nan(length(y),1);
% d = nan(length(kList),1);
% for kInd = 1:length(kList)
%     % Split train and test
%     te = k==kList(kInd);
% 
%     % Polar space normalization
%     x = polarSpaceNormalization(X,SVMspace,te);
% 
% %     % Within-session feature selection
% %     % get feature selection stats
% %     featStat = getFeatStat_ws(x,y,te,SVMspace);
% %     
% %     % apply feature selection
% %     switch SVMspace
% %         case {'cart_HT' 'cartNoAmp_HT' 'cartNoDelay_HT'...
% %                 'cartReal_T'...
% %                 'polMag_T'}
% %             featStat = abs(featStat);
% %             x = x(:,featStat>prctile(featStat,10));
% %         case {'cart' 'cartNoAmp' 'cartNoDelay'...
% %                 'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'...
% %                 'polMag'}
% %         otherwise
% %             error('X')
% %     end
% 
%     % Cocktail bank normalization
%     switch SVMspace
%         case {'hr' 'hrNoAmp'}
%             x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
%         case {'cart' 'cart_HT' 'cart_HTbSess'...
%                 'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'}
%             x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
%             x = cat(2,real(x),imag(x));
%         case {'cartReal' 'cartReal_T'}
%             x = real(x);
%             x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
%         case 'cartImag'
%             x = imag(x);
%             x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
%         case 'pol'
%             x = cat(2,angle(x),abs(x));
%             x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
%         case {'polMag' 'polMag_T' 'cartNoDelay' 'cartNoDelay_HT'  'cartNoDelay_HTbSess'}
%             x = abs(x);
%             x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
%         case 'polDelay'
%             x = angle(x);
%             x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
%         otherwise
%             error('x')
%     end
% 
%     % Train SVM
%     model = svmtrain(y(~te,:),x(~te,:),'-t 0 -q');
% %     w = model.sv_coef'*model.SVs;
% %     b = model.rho;
% %     yHat = cat(2,real(x(~te,:)),imag(x(~te,:)))*w';
%     
%     % Test SVM
% %     [yTr(~te,kInd), ~, yHatTr(~te,kInd)] = svmpredict(y(~te,:),x(~te,:),model,'-q');
%     [yTe(te,1), ~, yHatTe(te,1)] = svmpredict(y(te,:),x(te,:),model,'-q');
%     d(kInd) = size(x,2);
% end





function [x,nVox,nDim] = complex2svm(x,SVMspace)
% Output SVM ready data
nVox = size(x,2);
switch SVMspace
    case {'cart' 'cartNoAmp' 'cartNoAmp_affineRot' 'cartNoAmp_affineRot_affineCart' 'cart_roi' 'cart_affineRot'}
        x = cat(2,real(x),imag(x));
    case {'cartNoDelay' 'cartReal' 'cartReal_affineRot'}
        x = real(x);
    case {'cartNoAmpImag' 'cartImag' 'cartImag_affineRot' 'cartNoAmpImag_affineRot'}
        x = imag(x);
    otherwise
        error('X')
end
nDim = size(x,2);





% function featStat = getFeatStat_ws(x,y,te,SVMspace)
% nVox = size(x,2);
% switch SVMspace
%     case 'cart_HT'
%         x = [real(x) imag(x)];
%         
%         featStat = nan(1,nVox);
%         x = cat(1,x(~te & y==1,:),x(~te & y==2,:));
%         for voxInd = 1:nVox
%             stats = T2Hot2d(x(:,[0 nVox]+voxInd));
%             featStat(voxInd) = stats.T2;
%         end
%     case 'cartNoAmp_HT'
%         x = angle(x);
%         
%         featStat = nan(1,size(x,2));
%         for voxInd = 1:size(x,2)
%             [~, F] = circ_htest(x(~te & y==1,voxInd), x(~te & y==2,voxInd));
%             featStat(voxInd) = F;
%         end
%     case {'polMag_T' 'cartNoDelay_HT'}
%         x = abs(x);
%         
%         [~,~,~,STATS] = ttest(x(~te & y==1,:),x(~te & y==2,:));
%         featStat = STATS.tstat;
%     case 'cartReal_T'
%         x = real(x);
%         
%         [~,~,~,STATS] = ttest(x(~te & y==1,:),x(~te & y==2,:));
%         featStat = STATS.tstat;
%     case {'cart' 'cartNoAmp' 'cartNoDelay'...
%             'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'...
%             'polMag'}
%     otherwise
%         error('X')
% end



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


function res = perfMetric(y,yHat,k)
warning('off','stats:perfcurve:SubSampleWithMissingClasses')
averageWR = 1;
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
        'distT_p',[]);
    return
end

if averageWR
    nRun = length(unique(k))*length(unique(y));
    y = mean(reshape(y,[length(y)/nRun nRun]),1)';
    yHat = mean(reshape(yHat,[length(yHat)/nRun nRun]),1)';
end

res.y = {y};
res.yHat = {yHat};
res.nObs = length(y);
% acc
res.hit = sum((yHat<0)+1==y);
res.acc = res.hit./res.nObs;
[~,pci] = binofit(res.hit,res.nObs,0.1);
res.acc_CI5 = pci(1);
res.acc_CI95 = pci(2);
[~,pci] = binofit(res.nObs/2,res.nObs,0.1);
res.acc_thresh = pci(2);
res.acc_p = binocdf(res.hit,res.nObs,0.5,'upper');
% auc
[~,~,~,auc] = perfcurve(y,yHat,1,'NBOOT',2^10);
res.auc = auc(1);
res.auc_CI5 = auc(2);
res.auc_CI95 = auc(3);
% distT
[~,P,~,STATS] = ttest(yHat(y==1),yHat(y==2));
res.distT = STATS.tstat;
res.distT_p = P;


function [resSess,resSubj,resGroup] = summarizePerf(res_sess)
allField = fields(res_sess);
for i = 1:length(allField)
    if isnumeric(res_sess(1).(allField{i}))
        resSess.(allField{i}) = nan(size(res_sess));
        resSess.(allField{i})(:) = [res_sess.(allField{i})];
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
resSess = orderfields(resSess,[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 21 16 17 18 19 20]);
                              

resSubj.y = cell(size(resSess.y,1),1);
resSubj.yHat = cell(size(resSess.yHat,1),1);
for subjInd = 1:size(resSess.y,1)
    resSubj.y{subjInd} = cell2mat(resSess.y(subjInd,:)');
    resSubj.yHat{subjInd} = cell2mat(resSess.yHat(subjInd,:)');
end
resSubj.nObs = sum(resSess.nObs,2);
resSubj.hit = sum(resSess.hit,2);
resSubj.acc = resSubj.hit./resSubj.nObs;
[~,pci] = binofit(resSubj.hit,resSubj.nObs,0.1);
resSubj.acc_CI5 = pci(:,1);
resSubj.acc_CI95 = pci(:,2);
[~,pci] = binofit(resSubj.nObs/2,resSubj.nObs,0.1);
resSubj.acc_thresh = pci(:,2);
resSubj.acc_p = binocdf(resSubj.hit,resSubj.nObs,0.5,'upper');
resSubj.acc_fdr = mafdr(resSubj.acc_p,'BHFDR',true);
resSubj.auc = nan(size(resSubj.y));
resSubj.distT = nan(size(resSubj.y));
resSubj.distT_p = nan(size(resSubj.y));
for subjInd = 1:size(resSubj.y,1)
    [~,~,~,auc] = perfcurve(resSubj.y{subjInd},resSubj.yHat{subjInd},1);
    resSubj.auc(subjInd) = auc(1);
    resSubj.auc_CI5(subjInd) = nan;
    resSubj.auc_CI95(subjInd) = nan;
%     [~,~,~,auc] = perfcurve(resSubj.y{subjInd},resSubj.yHat{subjInd},1,'NBOOT',2^10);
%     resSubj.auc(subjInd) = auc(1);
%     resSubj.auc_CI5(subjInd) = auc(2);
%     resSubj.auc_CI95(subjInd) = auc(3);
    
    [~,P,~,STATS] = ttest(resSubj.yHat{subjInd}(resSubj.y{subjInd}==1),resSubj.yHat{subjInd}(resSubj.y{subjInd}==2));
    resSubj.distT(subjInd) = STATS.tstat;
    resSubj.distT_p(subjInd) = P;
end
resSubj.nVoxOrig = round(mean(resSess.nVoxOrig,2));
resSubj.nVox = round(mean(resSess.nVox,2));
resSubj.svmSpace = resSess.svmSpace(:,1);

resGroup.y = cell2mat(resSubj.y);
resGroup.yHat = cell2mat(resSubj.yHat);
resGroup.nObs = sum(resSubj.nObs,1);
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

[~,~,~,auc] = perfcurve(resGroup.y,resGroup.yHat,1,'NBOOT',2^11);
resGroup.auc = auc(1);
resGroup.auc_CI5 = auc(2);
resGroup.auc_CI95 = auc(3);
[~,P,~,STATS] = ttest(resSubj.auc,0.5,'tail','right');
resGroup.auc_T = STATS.tstat;
resGroup.auc_P = P;
[P,~,STATS] = signrank(resSubj.auc,0.5,'tail','right');
resGroup.auc_wilcoxonSignedrank = STATS.signedrank;
resGroup.auc_wilcoxonP = P;

[~,P,~,STATS] = ttest(resGroup.yHat(resGroup.y==1),resGroup.yHat(resGroup.y==2),'tail','right');
resGroup.distT = STATS.tstat;
resGroup.distT_p = P;
[~,P,~,STATS] = ttest(resSubj.distT,0,'tail','right');
resGroup.distT_T = STATS.tstat;
resGroup.distT_P = P;
resGroup.nVoxOrig = round(mean(resSubj.nVoxOrig,1));
resGroup.nVox = round(mean(resSubj.nVox,1));
resGroup.svmSpace = resSubj.svmSpace{1};

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