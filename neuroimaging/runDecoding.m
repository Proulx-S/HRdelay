function res = runDecoding(SVMspace,verbose,nPerm,figOption)
doAntiAntiLearning = 0;
if ~exist('verbose','var')
    verbose = 1;
end
if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end
if ~exist('SVMspace','var') || isempty(SVMspace)
    SVMspace = 'cart'; % 'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'
    % 'hr' 'hrNoAmp' 'cart' 'cartNoAmp' cartNoAmp_HT 'cartReal', 'cartImag', 'pol', 'polMag' 'polMag_T' or 'polDelay'
end
if isstruct(SVMspace)
    doPerm = 1;
    res = SVMspace;
    SVMspace = res.summary.SVMspace;
    if ~exist('nPerm','var') || isempty(nPerm)
        nPerm = 100;
    end
else
    doPerm = 0;
end

if doPerm
    filename = fullfile(pwd,mfilename);
    filename = fullfile(filename,[SVMspace '_' num2str(nPerm) 'perm']);
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

if ismac
    repoPath = '/Users/sebastienproulx/OneDrive - McGill University/dataBig';
else
    repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
end
dataDir = 'C-derived\DecodingHR';
funPath = fullfile(repoPath,dataDir,'fun');
funLevel = 'z';
fileSuffix = '_preprocAndShowMasks.mat';

%make sure everything is forward slash for mac, linux pc compatibility
for tmpPath = {'repoPath' 'dataDir' 'funPath'}
    eval([char(tmpPath) '(strfind(' char(tmpPath) ',''\''))=''/'';']);
end


%% Preload param
tmp = dir(fullfile(funPath,funLevel,'*.mat'));
for i = 1:length(tmp)
    curFile = fullfile(funPath,funLevel,tmp(i).name);
    load(curFile,'param')
    if exist('param','var')
        break
    end
end
if ~exist('param','var')
    error('Analysis parameters not found!')
end
subjList = param.subjList;


%% Load data
dCAll = cell(size(subjList,1),1);
for subjInd = 1:size(subjList,1)
    curFile = fullfile(funPath,funLevel,[subjList{subjInd} fileSuffix]);
    if verbose; disp(['loading: ' curFile]); end
    
    load(curFile,'dC');
    dCAll{subjInd} = dC;
end
dC = dCAll; clear dAll
sessList = fields(dC{1});

if verbose
    disp('---');
    disp(['SVM space: ' SVMspace]);
end

%% Reorganize
dP = cell(size(dC));
for subjInd = 1:length(dC)
    for sessInd = 1:length(sessList)
        sess = ['sess' num2str(sessInd)];
        dP{subjInd,sessInd} = dC{subjInd}.(sessList{sessInd});
        dC{subjInd}.(sessList{sessInd}) = [];
    end
end
clear dC

%% Between-session feature selection
featSel = cell(size(dP));
for i = 1:numel(dP)
    featSel{i}.ind = true(1,size(dP{i}.data,2));
    featSel{i}.info = 'V1';
    % Select non-vein voxels
    featSel{i}.ind = featSel{i}.ind & ~dP{i}.vein_mask;
    featSel{i}.info = strjoin({featSel{i}.info 'nonVein'},' & ');
    % Select active voxels
    featSel{i}.ind = featSel{i}.ind & dP{i}.anyCondActivation_mask;
    featSel{i}.info = strjoin({featSel{i}.info 'active'},' & ');
    % Select most discrimant voxels
    featSel{i}.ind = featSel{i}.ind & dP{i}.discrim_mask;
    featSel{i}.info = strjoin({featSel{i}.info 'mostDisciminant'},' & ');
end

% %% Between-session feature selection
% dP = cell(length(dC),2);
% featSel = repmat(struct('ind',[],'info',''),size(dP));
% for subjInd = 1:length(dC)
%     for sessInd = 1:length(sessList)
%         sess = ['sess' num2str(sessInd)];
%         %         ind.(sess) = true(1,size(dC{subjInd}.(sess).data,2));
%         featSel(subjInd,sessInd).ind = true(1,size(dC{subjInd}.(sess).data,2));
%         featSel(subjInd,sessInd).info = strjoin({featSel(subjInd,sessInd).info 'V1'},'');
%         % Select non-vein voxels
%         %         ind.(sess) = ind.(sess) & ~dC{subjInd}.(sess).vein_mask;
%         featSel(subjInd,sessInd).ind = featSel(subjInd,sessInd).ind & ~dC{subjInd}.(sess).vein_mask;
%         featSel(subjInd,sessInd).info = strjoin({featSel(subjInd,sessInd).info 'nonVein'},' & ');
%         % Select active voxels
%         %         ind.(sess) = ind.(sess) & dC{subjInd}.(sess).anyCondActivation_mask;
%         featSel(subjInd,sessInd).ind = featSel(subjInd,sessInd).ind & dC{subjInd}.(sess).anyCondActivation_mask;
%         featSel(subjInd,sessInd).info = strjoin({featSel(subjInd,sessInd).info 'active'},' & ');
%         % Select most discrimant voxels
%         %         ind.(sess) = ind.(sess) & dC{subjInd}.(sess).discrim_mask;
%         featSel(subjInd,sessInd).ind = featSel(subjInd,sessInd).ind & dC{subjInd}.(sess).discrim_mask;
%         featSel(subjInd,sessInd).info = strjoin({featSel(subjInd,sessInd).info 'mostDisciminant'},' & ');
%     end
%     
%     for sessInd = 1:length(sessList)
%         
%         % Apply voxel selection
%         for sessInd = 1:length(sessList)
%             sessFeat = ['sess' num2str(~(sessInd-1)+1)]; % the session on which voxel selection is defined
%             sess = ['sess' num2str(sessInd)]; % the session on which voxel selection is applied
%             
%             allFields = fields(dC{subjInd}.(sess));
%             nVox = size(dC{subjInd}.(sess).data,2);
%             for i = 1:length(allFields)
%                 if (isnumeric(dC{subjInd}.(sess).(allFields{i})) || islogical(dC{subjInd}.(sess).(allFields{i})) ) && size(dC{subjInd}.(sess).(allFields{i}),2)==nVox
%                     % remove unslected voxels identified from the other session
%                     dC{subjInd}.(sess).(allFields{i})(:,~ind.(sessFeat),:,:) = [];
%                 end
%             end
%             dP{subjInd,sessInd} = dC{subjInd}.(sessList{sessInd});
%         end
%         dC{subjInd} = [];
%     end
%     clear dC


%% Example plot of trigonometric (polar) representation
i = 1;
[~,b] = sort(dP{i}.anyCondActivation_F,'descend');
b = b(featSel{i}.ind);
f = plotPolNormExample(dP{i}.data(:,b(1),1:3),SVMspace);
if figOption.save
    filename = fullfile(pwd,mfilename);
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,SVMspace);
    f.Color = 'none';
    set(findobj(f.Children,'type','Axes'),'color','none')
    saveas(f,[filename '.svg']); if verbose; disp([filename '.svg']); end
    f.Color = 'w';
    set(findobj(f.Children,'type','Axes'),'color','w')
    saveas(f,[filename '.fig']); if verbose; disp([filename '.fig']); end
    saveas(f,[filename '.jpg']); if verbose; disp([filename '.jpg']); end
end

%% Some more independant-sample normalization
switch SVMspace
    case 'hrNoAmp' % scale each hr to 1 accroding to sin fit
        for i = 1:numel(dP)
            dP{i}.hr = dP{i}.hr./abs(dP{i}.data./100);
        end
    case {'cart' 'cart_HT' 'cart_HTbSess'...
            'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'...
            'cartNoDelay' 'cartNoDelay_HT' 'cartNoDelay_HTbSess'...
            'cartReal' 'cartReal_T'...
            'polMag' 'polMag_T'...
            'polDelay'}
        % Does not apply
    otherwise
        error('X')
end

%% Between-session svm
svmModel = cell(size(dP));
nrmlz = cell(size(dP));
yHat_tr = cell(size(dP));
svmModelK = cell(size(dP));
nrmlzK = cell(size(dP));
yHatK = cell(size(dP));
yHatK_tr = cell(size(dP));
yHat = cell(size(dP));
% polNorm = repmat(struct('rhoScale',[],'thetaShift',[]),size(dP));
% svmNorm = repmat(struct('scale',[],'shift',[]),size(dP));
% svmModel = repmat(struct('w',[],'b',[],'Paramters',[],'info',[]),size(dP));
% resTr.y  = cell(size(dP));
% resTr.yHat  = cell(size(dP));
% resTr.nObs = nan(size(dP));
% resTr.acc  = nan(size(dP));
% decision_valuesDiff__Tr = cell(size(dP));
% if doAntiAntiLearning
%     modelOneClass1 = cell(size(dP));
%     modelOneClass2 = cell(size(dP));
% end

% Within-session training and testing
Y = cell(size(dP));
for i = 1:numel(dP)
    [subjInd,sessInd] = ind2sub(size(dP),i);
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
        disp(['for sess ' num2str(i) ' of ' num2str(numel(dP))])
        tic
    end
    
    [X,y,k] = getXYK(dP{subjInd,sessInd},SVMspace);
    
    % SVM on within-session crossvalidation folds
    % feature selection
    x = X(:,featSel{subjInd,sessIndCross}.ind);
    % train
    [svmModelK{i},nrmlzK{i}] = SVMtrain(y,x,SVMspace,k);
    % test
    [yHatK{i},yHatK_tr{i}] = SVMtest(y,x,svmModelK{i},nrmlzK{i},k);
    yHatK{i} = nanmean(yHatK{i},2);
    
    
    % SVM on full sessions
    % feature-selected data
    x = X(:,featSel{subjInd,sessInd}.ind);
    % train
    [svmModel{i},nrmlz{i}] = SVMtrain(y,x,SVMspace);
    % test (not crossvalidated)
    [yHat_tr{i},~] = SVMtest(y,x,svmModel{i},nrmlz{i});
    
    
    % Output
    Y{i} = y;
end

% Between-session testing
for i = 1:numel(dP)
    [subjInd,testInd] = ind2sub(size(dP),i);
    switch testInd
        case 1
            trainInd = 2;
        case 2
            trainInd = 1;
        otherwise
            error('X')
    end
    % Get cross-session feature-selected data
    [X,y,~] = getXYK(dP{subjInd,testInd},SVMspace);
    x = X(:,featSel{subjInd,trainInd}.ind);
    
    % Test
    [yHat{subjInd,testInd},~] = SVMtest(y,x,svmModel{subjInd,trainInd});
end

%% Compute performance metrics
% Between-sess
resBS_sess = repmat(perfMetric,[size(dP,1) size(dP,2)]);
for i = 1:numel(dP)
    resBS_sess(i) = perfMetric(Y{i},yHat{i});
end
acc_fdr = mafdr([resBS_sess.acc_p],'BHFDR',true);
for i = 1:numel(dP)
    resBS_sess(i).acc_fdr = acc_fdr(i);
    resBS_sess(i).nVoxOrig = size(dP{i}.data,2);
    resBS_sess(i).nVox = sum(featSel{i}.ind);
    resBS_sess(i).SVMspace = SVMspace;
end
resBS_sess = orderfields(resBS_sess,[1 2 3 4 5 6 7 11 8 9 10 12 13 14]);
% Within-sess
resWS_sess = repmat(perfMetric,[size(dP,1) size(dP,2)]);
for i = 1:numel(dP)
    resWS_sess(i) = perfMetric(Y{i},yHatK{i});
end
acc_fdr = mafdr([resWS_sess.acc_p],'BHFDR',true);
for i = 1:numel(dP)
    resWS_sess(i).acc_fdr = acc_fdr(i);
    resWS_sess(i).nVoxOrig = size(dP{i}.data,2);
    resWS_sess(i).nVox = sum(featSel{i}.ind);
    resWS_sess(i).SVMspace = SVMspace;
end
resWS_sess = orderfields(resWS_sess,[1 2 3 4 5 6 7 11 8 9 10 12 13 14]);

%% Summarize group performances
[resBSsess,resBSsubj,resBSgroup] = summarizePerf(resBS_sess);
[resWSsess,resWSsubj,resWSgroup] = summarizePerf(resWS_sess);

figure('WindowStyle','docked');
metric = 'acc';
hPlot = plot(resWSsess.(metric)',resBSsess.(metric)','o'); hold on
for subjInd = 1:size(hPlot,1)
    hPlot(subjInd).MarkerFaceColor = hPlot(subjInd).Color;
    hPlot(subjInd).MarkerEdgeColor = 'k';
end
switch metric
    case {'acc' 'auc'}
        xlim([0 1])
        ylim([0 1])
        plot(xlim,[1 1].*0.5,'-k')
        plot([1 1].*0.5,ylim,'-k')
        ax = gca;
        ax.PlotBoxAspectRatio = [1 1 1];
    case 'distT'
        plot(xlim,[1 1].*0,'-k')
        plot([1 1].*0,ylim,'-k')
end
xlabel('Within-sess')
ylabel('Between-sess')

%% Print info
if verbose
    disp('-----------------')
    disp('*Between-session*')
    printRes(resBSsess,resBSsubj,resBSgroup)
    disp(' ')
    disp('*Within-session*')
    printRes(resWSsess,resWSsubj,resWSgroup)
    disp('-----------------')
end

% %% Add info
% if ~doPerm
%     
% else
%     res.perm.info = 'subj x sess x perm';
%     nObs = permute(repmat(res.nObs,[1 1 nPerm]),[3 1 2]);
%     res.perm.summary.hit = sum(res.perm.acc(:,:).*nObs(:,:),2);
%     res.perm.summary.nObs = sum(nObs(:,:),2);
%     res.perm.summary.acc = res.perm.summary.hit./res.perm.summary.nObs;
%     res.perm.summary.accThresh = prctile(res.perm.summary.acc,95);
%     res.perm.summary.p = sum(res.perm.summary.acc>res.summary.acc)./nPerm;
% 
%     res.perm.acc = permute(res.perm.acc,[2 3 1]);
%     res.perm.auc = permute(res.perm.auc,[2 3 1]);
%     res.perm.distT = permute(res.perm.distT,[2 3 1]);
% 
%     if verbose
%         disp('Group results:')
%         disp(['  hit    =' num2str(res.summary.hit) '/' num2str(res.summary.nObs)])
%         disp(['  acc    =' num2str(res.summary.acc*100,'%0.2f%%')])
%         disp(' permutation test stats')
%         disp(['  thresh =' num2str(res.perm.summary.accThresh*100,'%0.2f%%')])
%         disp(['  p      =' num2str(res.perm.summary.p,'%0.3f')])
%     end
% 
%     filename = fullfile(pwd,mfilename);
%     if ~exist(filename,'dir'); mkdir(filename); end
%     filename = fullfile(filename,[SVMspace '_' num2str(nPerm) 'perm']);
%     save(filename,'res')
%     if verbose; disp([filename '.mat']); end
% end

function [x,y,k] = getXYK(dP,SVMspace)
% Define x(data), y(label) and k(xValFolds)
switch SVMspace
    case {'hr' 'hrNoAmp'}
        error('double-check that')
        nSamplePaire = size(dP.hr,1);
        x1 = dP.hr(:,:,1,:); x1 = x1(:,:);
        x2 = dP.hr(:,:,2,:); x2 = x2(:,:);
    case {'cart' 'cart_HT' 'cart_HTbSess'...
            'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'...
            'cartNoDelay' 'cartNoDelay_HT' 'cartNoDelay_HTbSess'...
            'cartReal' 'cartReal_T'...
            'polMag' 'polMag_T'...
            'polDelay'}
        nSamplePaire = size(dP.data,1);
        x1 = dP.data(:,:,1);
        x2 = dP.data(:,:,2);
    otherwise
        error('X')
end
y1 = 1.*ones(nSamplePaire,1);
k1 = (1:nSamplePaire)';
y2 = 2.*ones(nSamplePaire,1);
k2 = (1:nSamplePaire)';

x = cat(1,x1,x2); clear x1 x2
y = cat(1,y1,y2); clear y1 y2
k = cat(1,k1,k2); clear k1 k2

function [svmModel,normTr] = SVMtrain(y,x,SVMspace,k)
Y = y;
X = x;
if ~exist('k','var') || isempty(k) || all(k(:)==0)
    k = ones(size(y)); % no crossvalidation
    kList = 0;
else
    kList = unique(k);
end


normTr.k = [];
normTr.svmSpace = '';
normTr.pol = struct('rhoScale',[],'thetaShift',[]);
normTr.svm = struct('scale',[],'shift',[]);
normTr = repmat(normTr,[1 1 length(kList)]);
svmModel = repmat(struct('k',[],'svmSpace',[],'w',[],'b',[],'Paramters',[],'info',[]),[1 1 length(kList)]);
for kInd = 1:length(kList)
    % Split train and test
    te = k==kList(kInd);
    y = Y(~te,:);
    x = X(~te,:);
    
    % Normalize
    [x,normTr(:,:,kInd).pol] = polarSpaceNormalization(x,SVMspace);
    [x,normTr(:,:,kInd).svm] = svmSpaceNormalization(x,SVMspace);
    normTr(:,:,kInd).k = kList(kInd);
    normTr(:,:,kInd).svmSpace = SVMspace;
    
    % Train
    model = svmtrain(y,x,'-t 0 -q');
    w = model.sv_coef'*model.SVs;
    b = model.rho;
    
    % Store
    svmModel(:,:,kInd).k = kList(kInd);
    svmModel(:,:,kInd).w = w;
    svmModel(:,:,kInd).b = b;
    svmModel(:,:,kInd).Paramters = model.Parameters;
    svmModel(:,:,kInd).info = '   yHat = x*w''-b   ';
    svmModel(:,:,kInd).svmSpace = SVMspace;
end

function [yHat,yHatTr] = SVMtest(y,x,svmModel,nrmlz,k)
X = x;
if ~exist('k','var') || isempty(k) || all(k(:)==0)
    k = ones(size(y)); % no crossvalidation
end
kList = [svmModel.k];

yHat = nan([size(y,1) length(kList)]);
yHatTr = nan([size(y,1) length(kList)]);
for kInd = 1:length(kList)
    x = X;
    % Split train and test
    if kList(kInd)==0
        te = true(size(y));
    else
        te = k==kList(kInd);
    end
    
    % Normalize
    if ~exist('nrmlz','var') || isempty(nrmlz)
        [x,~] = polarSpaceNormalization(x,svmModel(kInd).svmSpace);
        [x,~] = svmSpaceNormalization(x,svmModel(kInd).svmSpace);
    else
        [x,~] = polarSpaceNormalization(x,nrmlz(kInd),te);
        [x,~] = svmSpaceNormalization(x,nrmlz(kInd),te);
    end
    
    w = svmModel(kInd).w;
    b = svmModel(kInd).b;
    yHat(te,kInd) = x(te,:)*w'-b;
    yHatTr(~te,kInd) = x(~te,:)*w'-b;
end
    


function [yTe,d,yHatTe] = xValSVM(x,y,k,SVMspace)
X = x;
kList = unique(k);
% yTr = nan(length(y),length(kList));
yTe = nan(length(y),1);
% yHatTr = nan(length(y),length(kList));
yHatTe = nan(length(y),1);
d = nan(length(kList),1);
for kInd = 1:length(kList)
    % Split train and test
    te = k==kList(kInd);

    % Polar space normalization
    x = polarSpaceNormalization(X,SVMspace,te);

%     % Within-session feature selection
%     % get feature selection stats
%     featStat = getFeatStat_ws(x,y,te,SVMspace);
%     
%     % apply feature selection
%     switch SVMspace
%         case {'cart_HT' 'cartNoAmp_HT' 'cartNoDelay_HT'...
%                 'cartReal_T'...
%                 'polMag_T'}
%             featStat = abs(featStat);
%             x = x(:,featStat>prctile(featStat,10));
%         case {'cart' 'cartNoAmp' 'cartNoDelay'...
%                 'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'...
%                 'polMag'}
%         otherwise
%             error('X')
%     end

    % Cocktail bank normalization
    switch SVMspace
        case {'hr' 'hrNoAmp'}
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
        case {'cart' 'cart_HT' 'cart_HTbSess'...
                'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'}
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
            x = cat(2,real(x),imag(x));
        case {'cartReal' 'cartReal_T'}
            x = real(x);
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
        case 'cartImag'
            x = imag(x);
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
        case 'pol'
            x = cat(2,angle(x),abs(x));
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
        case {'polMag' 'polMag_T' 'cartNoDelay' 'cartNoDelay_HT'  'cartNoDelay_HTbSess'}
            x = abs(x);
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
        case 'polDelay'
            x = angle(x);
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
        otherwise
            error('x')
    end

    % Train SVM
    model = svmtrain(y(~te,:),x(~te,:),'-t 0 -q');
%     w = model.sv_coef'*model.SVs;
%     b = model.rho;
%     yHat = cat(2,real(x(~te,:)),imag(x(~te,:)))*w';
    
    % Test SVM
%     [yTr(~te,kInd), ~, yHatTr(~te,kInd)] = svmpredict(y(~te,:),x(~te,:),model,'-q');
    [yTe(te,1), ~, yHatTe(te,1)] = svmpredict(y(te,:),x(te,:),model,'-q');
    d(kInd) = size(x,2);
end

function [x,polNorm] = polarSpaceNormalization(x,SVMspace,te)
if ~exist('te','var') || isempty(te)
    te = false(size(x,1),1);
end
if isstruct(SVMspace)
    polNorm = SVMspace.pol;
    SVMspace = SVMspace.svmSpace;
    computeNorm = false;
else
    computeNorm = true;
end
sz = size(x);
nDim = ndims(x);
switch nDim
    case 2
    case 3
        % Redim
        xTmp = nan(prod(sz([1 3])),sz(2));
        teTmp = nan(prod(sz([1 3])),1);
        for i = 1:sz(3)
            xTmp( (1:sz(1)) + sz(1)*(i-1) , : ) = x(:,:,i);
            teTmp( (1:sz(1)) + sz(1)*(i-1) , 1 ) = false(sz(1),1);
        end
        x = xTmp;
        te = teTmp;
    otherwise
        error('dim1 must be samples and dim2 voxels, that''s it')
end

% Compute norm
switch SVMspace
    case {'hr' 'hrNoAmp'}
        % does not apply
    case {'cart' 'cart_HT' 'cart_HTbSess'...
            'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'...
            'cartNoDelay' 'cartNoDelay_HT' 'cartNoDelay_HTbSess'...
            'cartReal' 'cartReal_T'...
            'polMag' 'polMag_T'...
            'polDelay'}
        % rho
        switch SVMspace
            case {'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'...
                    'polDelay'}
                if computeNorm
                    polNorm.rhoScale = ones([1 size(x,2)]);
                end
                rho = ones(size(x));
            case {'cart' 'cart_HT' 'cart_HTbSess'...
                    'cartNoDelay' 'cartNoDelay_HT' 'cartNoDelay_HTbSess'...
                    'cartReal' 'cartReal_T'...
                    'polMag' 'polMag_T'}
                if computeNorm
                   polNorm.rhoScale = abs(mean(mean(x(~te,:),1),2)); 
                end
                rho = abs(x)./polNorm.rhoScale;
            otherwise
                error('X')
        end
        % thetaShift
        switch SVMspace
            case {'cart' 'cart_HT' 'cart_HTbSess'...
                    'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'...
                    'cartReal' 'cartReal_T'...
                    'polDelay'}
                if computeNorm
                    polNorm.thetaShift = angle(mean(mean(x(~te,:),1),2));
                end
                theta = angle(x) - polNorm.thetaShift; theta = wrapToPi(theta);
            case {'cartNoDelay' 'cartNoDelay_HT' 'cartNoDelay_HTbSess'...
                    'polMag' 'polMag_T'}
                if computeNorm
                    polNorm.thetaShift = zeros([1,size(x,2)]);
                end
                theta = zeros(size(x));
            otherwise
                error('X')
        end
        [u,v] = pol2cart(theta,rho); clear theta rho
        x = complex(u,v); clear u v
    otherwise
        error('X')
end

if nDim==3
    % Redim
    xTmp = nan(sz);
    for i = 1:sz(3)
        xTmp(:,:,i) = x( (1:sz(1)) + sz(1)*(i-1) , : );
    end
    x = xTmp;
end

function [x,svmNorm] = svmSpaceNormalization(x,SVMspace,te)
if ~exist('te','var') || isempty(te)
    te = false(size(x,1),1);
end
if isstruct(SVMspace)
    svmNorm = SVMspace.svm;
    SVMspace = SVMspace.svmSpace;
    computeNorm = false;
else
    computeNorm = true;
end

% Move to SVMspace and Cocktail bank normalize
switch SVMspace
    case {'hr' 'hrNoAmp'}
        error('code that')
        x = x-mean(x(~te,:),1);
    case {'cart' 'cart_HT' 'cart_HTbSess'...
            'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'}
        if computeNorm
            svmNorm.scale = ones(1,size(x(~te,:),2));
            svmNorm.shift = mean(x(~te,:),1);
        end
        x = (x-svmNorm.shift) ./ svmNorm.scale;
        x = cat(2,real(x),imag(x));
    case {'cartReal' 'cartReal_T'}
        x = real(x);
        if computeNorm
            svmNorm.scale = ones(1,size(x(~te,:),2));
            svmNorm.shift = mean(x(~te,:),1);
        end
        x = (x-svmNorm.shift) ./ svmNorm.scale;
    case 'cartImag'
        error('code that')
        x = imag(x);
        if computeNorm
            svmNorm.scale = ones(1,size(x(~te,:),2));
            svmNorm.shift = mean(x(~te,:),1);
        end
        x = (x-svmNorm.shift) ./ svmNorm.scale;
    case 'pol'
        error('code that')
        x = cat(2,angle(x),abs(x));
        x = x-mean(x(~te,:),1);
    case {'polMag' 'polMag_T' 'cartNoDelay' 'cartNoDelay_HT'  'cartNoDelay_HTbSess'}
        x = abs(x);
        if computeNorm
            svmNorm.scale = ones(1,size(x(~te,:),2));
            svmNorm.shift = mean(x(~te,:),1);
        end
        x = (x-svmNorm.shift) ./ svmNorm.scale;
    case 'polDelay'
        if computeNorm
            svmNorm.scale = ones(1,size(x(~te,:),2));
            svmNorm.shift = zeros(1,size(x(~te,:),2));
        end
        x = angle(x);
        x = (x - svmNorm.shift) ./ svmNorm.scale;
    otherwise
        error('x')
end


function featStat = getFeatStat_ws(x,y,te,SVMspace)
nVox = size(x,2);
switch SVMspace
    case 'cart_HT'
        x = [real(x) imag(x)];
        
        featStat = nan(1,nVox);
        x = cat(1,x(~te & y==1,:),x(~te & y==2,:));
        for voxInd = 1:nVox
            stats = T2Hot2d(x(:,[0 nVox]+voxInd));
            featStat(voxInd) = stats.T2;
        end
    case 'cartNoAmp_HT'
        x = angle(x);
        
        featStat = nan(1,size(x,2));
        for voxInd = 1:size(x,2)
            [~, F] = circ_htest(x(~te & y==1,voxInd), x(~te & y==2,voxInd));
            featStat(voxInd) = F;
        end
    case {'polMag_T' 'cartNoDelay_HT'}
        x = abs(x);
        
        [~,~,~,STATS] = ttest(x(~te & y==1,:),x(~te & y==2,:));
        featStat = STATS.tstat;
    case 'cartReal_T'
        x = real(x);
        
        [~,~,~,STATS] = ttest(x(~te & y==1,:),x(~te & y==2,:));
        featStat = STATS.tstat;
    case {'cart' 'cartNoAmp' 'cartNoDelay'...
            'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'...
            'polMag'}
    otherwise
        error('X')
end

function f = plotPolNormExample(x,SVMspace)
f = figure('WindowStyle','docked');
subplot(1,2,1); clear hPP
% polarplot(angle(x(:)),abs(x(:)),'.'); hold on
for condInd = 1:3
    hPP(condInd) = polarplot(angle(x(:,:,condInd)),abs(x(:,:,condInd)),'o'); hold on
    hPP(condInd).MarkerFaceColor = hPP(condInd).Color;
    hPP(condInd).MarkerEdgeColor = 'w';
    hPP(condInd).MarkerSize = 3.5;
end
drawnow
hPP1 = hPP; clear hPP
ax1 = gca;

% Normalize
x = polarSpaceNormalization(x,SVMspace);

% Plot after
subplot(1,2,2);
% polarplot(angle(xAfter(:)),abs(xAfter(:)),'.'); hold on
for condInd = 1:3
    hPP(condInd) = polarplot(angle(x(:,:,condInd)),abs(x(:,:,condInd)),'o'); hold on
    hPP(condInd).MarkerFaceColor = hPP(condInd).Color;
    hPP(condInd).MarkerEdgeColor = 'w';
    hPP(condInd).MarkerSize = 3.5;
end
drawnow
hPP2 = hPP; clear hPP
ax2 = gca;

ax = ax1;
ax.ThetaTickLabel = 12-ax.ThetaTick(1:end)/360*12;
ax.ThetaTickLabel(1,:) = '0 ';
ax.ThetaAxis.Label.String = {'delay' '(sec)'};
ax.ThetaAxis.Label.Rotation = 0;
ax.ThetaAxis.Label.HorizontalAlignment = 'left';
ax.RAxis.Label.String = 'amp (%BOLD)';
ax.RAxis.Label.Rotation = 80;
ax.Title.String = 'before';

ax = ax2;
ax.ThetaTickLabel = (-wrapTo180(ax.ThetaTick(1:end))/360*12);
ax.ThetaTickLabel(1,:) = '0 ';
% ax.ThetaAxis.Label.String = {'delay' '(sec)'};
% ax.ThetaAxis.Label.Rotation = 0;
% ax.ThetaAxis.Label.HorizontalAlignment = 'left';
ax.Title.String = 'after';

hSup = suptitle({'Polar space normalization' SVMspace});
hSup.Interpreter = 'none';
drawnow

function res = perfMetric(y,yHat)
if ~exist('y','var')
    res = struct(...
        'y',[],...
        'yHat',[],...
        'nObs',[],...
        'hit',[],...
        'acc',[],...
        'acc_thresh',[],...
        'acc_p',[],...);
        'auc',[],...
        'distT',[],...
        'distT_p',[]);
    return
end
res.y = {y};
res.yHat = {yHat};
res.nObs = length(y);
% acc
res.hit = sum((yHat<0)+1==y);
res.acc = res.hit./res.nObs;
[~,pci] = binofit(res.nObs/2,res.nObs,0.1);
res.acc_thresh = pci(2);
res.acc_p = binocdf(res.hit,res.nObs,0.5,'upper');
% auc
[~,~,~,res.auc] = perfcurve(y,yHat,1);
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
    else
        error('code that')
    end
end
resSess.info = 'subj x sess';
resSess.distT_fdr = resSess.distT_p;
resSess.distT_fdr(:) = mafdr(resSess.distT_p(:),'BHFDR',true);
resSess = orderfields(resSess,[1 2 3 4 5 6 7 8 9 10 11 16 12 13 14 15]);

resSubj.y = cell(size(resSess.y,1),1);
resSubj.yHat = cell(size(resSess.yHat,1),1);
for subjInd = 1:size(resSess.y,1)
    resSubj.y{subjInd} = cell2mat(resSess.y(subjInd,:)');
    resSubj.yHat{subjInd} = cell2mat(resSess.yHat(subjInd,:)');
end
resSubj.nObs = sum(resSess.nObs,2);
resSubj.hit = sum(resSess.hit,2);
resSubj.acc = resSubj.hit./resSubj.nObs;
[~,pci] = binofit(resSubj.nObs/2,resSubj.nObs,0.1);
resSubj.acc_thresh = pci(:,2);
resSubj.acc_p = binocdf(resSubj.hit,resSubj.nObs,0.5,'upper');
resSubj.acc_fdr = mafdr(resSubj.acc_p,'BHFDR',true);
resSubj.auc = nan(size(resSubj.y));
resSubj.distT = nan(size(resSubj.y));
resSubj.distT_p = nan(size(resSubj.y));
for subjInd = 1:size(resSubj.y,1)
    [~,~,~,resSubj.auc(subjInd)] = perfcurve(resSubj.y{subjInd},resSubj.yHat{subjInd},1);
    [~,P,~,STATS] = ttest(resSubj.yHat{subjInd}(resSubj.y{subjInd}==1),resSubj.yHat{subjInd}(resSubj.y{subjInd}==2));
    resSubj.distT(subjInd) = STATS.tstat;
    resSubj.distT_p(subjInd) = P;
end
resSubj.nVoxOrig = round(mean(resSess.nVoxOrig,2));
resSubj.nVox = round(mean(resSess.nVox,2));
resSubj.SVMspace = resSess.SVMspace(:,1);

resGroup.y = cell2mat(resSubj.y);
resGroup.yHat = cell2mat(resSubj.yHat);
resGroup.nObs = sum(resSubj.nObs,1);
resGroup.hit = sum(resSubj.hit,1);

resGroup.acc = resGroup.hit./resGroup.nObs;
[~,pci] = binofit(resGroup.nObs/2,resGroup.nObs,0.1);
resGroup.acc_thresh = pci(:,2);
resGroup.acc_p = binocdf(resGroup.hit,resGroup.nObs,0.5,'upper');
[~,P,~,STATS] = ttest(resSubj.acc,0.5,'tail','right');
resGroup.acc_T = STATS.tstat;
resGroup.acc_P = P;
[P,~,STATS] = signrank(resSubj.acc,0.5,'tail','right');
resGroup.acc_wilcoxonSignedrank = STATS.signedrank;
resGroup.acc_wilcoxonP = P;

[~,~,~,auc] = perfcurve(resGroup.y,resGroup.yHat,1,'NBOOT',1000);
resGroup.auc = auc(1);
resGroup.auc_CI = auc(2:3);
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
resGroup.SVMspace = resSubj.SVMspace{1};

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
disp(['  auc    =' num2str(resGroup.auc*100,'%0.2f') '; 90%CI=' num2str(resGroup.auc_CI,'%0.3f ')])
disp(['  distT  =' num2str(resGroup.distT,'%0.2f') '; p=' num2str(resGroup.distT_p,'%0.3f ')])
disp([' -randomEffect'])
disp(['  acc=' num2str(mean(resSubj.acc)*100,'%0.2f%%')])
disp(['  -student'])
disp(['   T=' num2str(resGroup.acc_T,'%0.2f') '; P=' num2str(resGroup.acc_P,'%0.3f')])
disp(['  -wilcoxon'])
disp(['   sRank=' num2str(resGroup.acc_wilcoxonSignedrank,'%0.2f') '; P=' num2str(resGroup.acc_wilcoxonP,'%0.3f')])

