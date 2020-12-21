function res = runDecoding(SVMspace,nPerm,figOption)
if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end
if ~exist('SVMspace','var') || isempty(SVMspace)
    SVMspace = 'cartNoAmp_HTbSess'; % 'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'
    % 'hr' 'hrNoAmp' 'cart' 'cartNoAmp' cartNoAmp_HT 'cartReal', 'cartImag', 'pol', 'polMag' 'polMag_T' or 'polDelay'

    % cartNoDelay_HT WIERD!!!
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

disp('---');
disp(['SVM space: ' SVMspace]);
if doPerm
    filename = fullfile(pwd,mfilename);
    filename = fullfile(filename,[SVMspace '_' num2str(nPerm) 'perm']);
    if exist([filename '.mat'],'file')
        load(filename)
        disp('permutation found on disk, skipping')
        disp('Group results:')
        disp(['  hit    =' num2str(res.summary.hit) '/' num2str(res.summary.nObs)])
        disp(['  acc    =' num2str(res.summary.acc*100,'%0.2f%%')])
        disp(' permutation test stats')
        disp(['  thresh =' num2str(res.perm.summary.accThresh*100,'%0.2f%%')])
        disp(['  p      =' num2str(res.perm.summary.p,'%0.3f') ')'])
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
    disp(['loading: ' curFile]);
    
    load(curFile,'dC');
    dCAll{subjInd} = dC;
end
dC = dCAll; clear dAll
sessList = fields(dC{1});


%% Between-session feature selection
dP = cell(length(dC),2);
for subjInd = 1:length(dC)
    for sessInd = 1:length(sessList)
        sess = ['sess' num2str(sessInd)];
        ind.(sess) = true(1,size(dC{subjInd}.(sess).data,2));
        % Select non-vein voxels
        ind.(sess) = ind.(sess) & ~dC{subjInd}.(sess).vein_mask;
        % Select active voxels
        ind.(sess) = ind.(sess) & dC{subjInd}.(sess).anyCondActivation_mask;
        % Select most discrimant voxels
        switch SVMspace
            case {'cart' 'cart_HT' 'cartNoAmp' 'cartNoAmp_HT' 'cartNoDelay' 'cartNoDelay_HT' 'cartReal_T' 'polMag' 'polMag_T'}
            case {'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'}
                x = dC{subjInd}.(sess).data(:,:,1:2);
                switch SVMspace
                    case 'cart_HTbSess'
                    case 'cartNoAmp_HTbSess'
                        [u,v] = pol2cart(angle(x),1);
                        x = complex(u,v);
                    case 'cartNoDelay_HTbSess'
                        [u,v] = pol2cart(0,abs(x));
                        x = complex(u,v);
                    otherwise
                        error('X')
                end
                y = []; for condInd = 1:size(x,3); y = cat(3,y,ones(size(x,1),1).*condInd); end
                x2 = []; y2 = [];
                for condInd = 1:size(x,3)
                    x2 = cat(1,x2,x(:,:,condInd));
                    y2 = cat(1,y2,y(:,:,condInd));
                end
                x = x2; y = y2; clear x2 y2
                switch SVMspace
                    case {'cart_HTbSess' 'cartNoAmp_HTbSess'}
                        % hotelling's Tsquared for cond1 vs cond2 on [sin cos]
                        featStat = nan(1,size(x,2));
                        for voxInd = 1:size(x,2)
                            xTmp = cat(1,x(y==1,voxInd),x(y==2,voxInd));
                            stats = T2Hot2d([real(xTmp) imag(xTmp)]);
                            featStat(voxInd) = stats.T2;
                        end
                    case 'cartNoDelay_HTbSess'
                        % ttest for cond1 vs cond2 on sin (real)
                        [~,~,~,STATS] = ttest(real(x(y==1,:)),real(x(y==2,:)));
                        featStat = STATS.tstat;
                    otherwise
                        error('X')
                end
                % select top 90% most discriminant (active and non-vein) voxels
                ind.(sess) = ind.(sess) & featStat>prctile(featStat(ind.(sess)),10);
            otherwise
                error('X')
        end
    end
    % Apply voxel selection
    for sessInd = 1:length(sessList)
        sessFeat = ['sess' num2str(~(sessInd-1)+1)]; % the session on which voxel selection is defined
        sess = ['sess' num2str(sessInd)]; % the session on which voxel selection is applied

        allFields = fields(dC{subjInd}.(sess));
        nVox = size(dC{subjInd}.(sess).data,2);
        for i = 1:length(allFields)
            if (isnumeric(dC{subjInd}.(sess).(allFields{i})) || islogical(dC{subjInd}.(sess).(allFields{i})) ) && size(dC{subjInd}.(sess).(allFields{i}),2)==nVox
                % remove unslected voxels identified from the other session
                dC{subjInd}.(sess).(allFields{i})(:,~ind.(sessFeat),:,:) = [];
            end
        end
        dP{subjInd,sessInd} = dC{subjInd}.(sessList{sessInd});
    end
    dC{subjInd} = [];
end
clear dC


%% Example plot of trigonometric (polar) representation
[~,b] = max(dP{1}.anyCondActivation_F);
f = plotPolNormExample(dP{1}.data(:,b,1:3),SVMspace);
if figOption.save
    filename = fullfile(pwd,mfilename);
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,SVMspace);
    f.Color = 'none';
    set(findobj(f.Children,'type','Axes'),'color','none')
    saveas(f,[filename '.svg']); disp([filename '.svg'])
    f.Color = 'w';
    set(findobj(f.Children,'type','Axes'),'color','w')
    saveas(f,[filename '.fig']); disp([filename '.fig'])
    saveas(f,[filename '.jpg']); disp([filename '.jpg'])
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
            'cartReal_T'...
            'polMag' 'polMag_T'}
        % Does not apply
    otherwise
        error('X')
end


%% Run svm
if ~doPerm
    res.nVox = nan(size(dP));
    res.nDim = nan(size(dP));
    res.acc  = nan(size(dP));
    res.auc  = nan(size(dP));
    res.distT = nan(size(dP));
    res.nObs = nan(size(dP));
    res.p = nan(size(dP));
    res.subjList = subjList;
else
    res.perm.acc = nan([nPerm size(res.acc)]);
    res.perm.auc = nan([nPerm size(res.auc)]);
    res.perm.distT = nan([nPerm size(res.distT)]);
end
for i = 1:numel(dP)
    if doPerm
        disp(['for sess ' num2str(i) ' of ' num2str(numel(dP))])
        tic
    end
    % Define x(data), y(label) and k(xValFolds)
    switch SVMspace
        case {'hr' 'hrNoAmp'}
            nSamplePaire = size(dP{i}.hr,1);
            x1 = dP{i}.hr(:,:,1,:); x1 = x1(:,:);
            x2 = dP{i}.hr(:,:,2,:); x2 = x2(:,:);
        case {'cart' 'cart_HT' 'cart_HTbSess'...
                'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'...
                'cartNoDelay' 'cartNoDelay_HT' 'cartNoDelay_HTbSess'...
                'cartReal_T'...
                'polMag' 'polMag_T'}
            nSamplePaire = size(dP{i}.data,1);
            x1 = dP{i}.data(:,:,1);
            x2 = dP{i}.data(:,:,2);
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

    % SVM
    if ~doPerm
        %cross-validated SVM
        [yTe,d,yHatTe] = xValSVM(x,y,k,SVMspace);
    else
        % with permutations
        yTe = nan(length(y),nPerm);
        yHatTe = nan(length(y),nPerm);
        y1 = y(y==1);
        y2 = y(y==2);
        parfor permInd = 1:nPerm
            %permute labels
            y = cat(2,y1,y2);
            for sInd = 1:size(y,1); y(sInd,:) = y(sInd,randperm(2)); end
            y = cat(1,y(:,1),y(:,2));
            %cross-validated SVM
            [yTe(:,permInd),~,yHatTe(:,permInd)] = xValSVM(x,y,k,SVMspace);
        end
        y = cat(1,y1,y2); clear y1 y2
    end

    % Evaluate SVM
    if ~doPerm
        res.nObs(i) = length(y);
        res.acc(i) = sum(yTe==y)./res.nObs(i);
        [FP,TP,T,AUC] = perfcurve(y,yHatTe,1);
        res.auc(i) = AUC;
%         figure('WindowStyle','docked')
%         plot(FP.*res.nObs(i)./2,TP.*res.nObs(i)./2,'r','linewidth',2)
%         xlabel('False Positives')
%         ylabel('True Positives')
%         ax = gca; ax.PlotBoxAspectRatio = [1 1 1];
        [~,~,~,STATS] = ttest(yHatTe(y==1),yHatTe(y==2));
        res.distT(i) = STATS.tstat;
        res.p(i) = binocdf(res.acc(i).*res.nObs(i),res.nObs(i),0.5,'upper');
        res.nDim(i) = mean(d);
        switch SVMspace
            case 'hr'
                res.nVox1(i) = size(dP{i}.hr,2); % before within-session feature selection
                res.nVox2(i) = mean(d)./size(dP{i}.hr,4); % before within-session feature selection
            case {'cart' 'cart_HT' 'cart_HTbSess'...
                    'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'...
                    'cartNoDelay' 'cartNoDelay_HT' 'cartNoDelay_HTbSess'}
                res.nVox1(i) = size(dP{i}.data,2); % before within-session feature selection
                res.nVox2(i) = mean(d)./2; % before within-session feature selection
            case {'polMag' 'polMag_T' 'cartReal_T'}
                res.nVox1(i) = size(dP{i}.data,2); % before within-session feature selection
                res.nVox2(i) = mean(d); % before within-session feature selection
            otherwise
                error('X')
        end
    else
        res.perm.acc(:,i) = sum(yTe==y,1)./res.nObs(i);
        for permInd = 1:nPerm
            [~,~,~,res.auc(permInd,i)] = perfcurve(y,yHatTe(:,permInd),1);
            [~,~,~,STATS] = ttest(yHatTe(y==1,permInd),yHatTe(y==2,permInd));
            res.distT(permInd,i) = STATS.tstat;
        end
    end
    if doPerm
        toc
    end
end
% figure('WindowStyle','docked')
% x = res.auc(:);
% y = res.acc(:);
% scatter(x,y); hold on
% ax = gca;
% ax.PlotBoxAspectRatio = [1 1 1];
% grid on
% xlim([0 1]); ylim([0 1]);
% uistack(plot([0 1],[0 1],'k'),'bottom')


%% Add info
if ~doPerm
    res.info = 'subj x sess';
    res.summary.SVMspace = SVMspace;
    res.summary.hit = sum(res.acc(:).*res.nObs(:));
    res.summary.nObs = sum(res.nObs(:));
    res.summary.acc = res.summary.hit/res.summary.nObs;
    res.summary.auc = mean(mean(res.auc,2),1);
    res.summary.distT = mean(mean(res.distT,2),1);
    [~,pci] = binofit(res.summary.nObs/2,res.summary.nObs,0.1);
    res.summary.accThresh = pci(2);
    res.summary.p = binocdf(res.summary.hit,res.summary.nObs,0.5,'upper');
    disp('Group results:')
    disp(['  hit    =' num2str(res.summary.hit) '/' num2str(res.summary.nObs)])
    disp(['  acc    =' num2str(res.summary.acc*100,'%0.2f%%')])
    disp(['  auc    =' num2str(res.summary.auc*100,'%0.2f')])
    disp(['  distT  =' num2str(res.summary.distT,'%0.2f')])
    disp(' binomial stats')
    disp(['  thresh =' num2str(res.summary.accThresh*100,'%0.2f')])
    disp(['  p      =' num2str(res.summary.p,'%0.3f') ')'])
else
    res.perm.info = 'subj x sess x perm';
    nObs = permute(repmat(res.nObs,[1 1 nPerm]),[3 1 2]);
    res.perm.summary.hit = sum(res.perm.acc(:,:).*nObs(:,:),2);
    res.perm.summary.nObs = sum(nObs(:,:),2);
    res.perm.summary.acc = res.perm.summary.hit./res.perm.summary.nObs;
    res.perm.summary.accThresh = prctile(res.perm.summary.acc,95);
    res.perm.summary.p = sum(res.perm.summary.acc>res.summary.acc)./nPerm;

    res.perm.acc = permute(res.perm.acc,[2 3 1]);
    res.perm.auc = permute(res.perm.auc,[2 3 1]);
    res.perm.distT = permute(res.perm.distT,[2 3 1]);

    disp('Group results:')
    disp(['  hit    =' num2str(res.summary.hit) '/' num2str(res.summary.nObs)])
    disp(['  acc    =' num2str(res.summary.acc*100,'%0.2f%%')])
    disp(' permutation test stats')
    disp(['  thresh =' num2str(res.perm.summary.accThresh*100,'%0.2f%%')])
    disp(['  p      =' num2str(res.perm.summary.p,'%0.3f') ')'])

    filename = fullfile(pwd,mfilename);
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,[SVMspace '_' num2str(nPerm) 'perm']);
    save(filename,'res')
    disp([filename '.mat'])
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
    x = polarSpaceNormalization(X,te,SVMspace);

    % Within-session feature selection
    % get feature selection stats
    featStat = getFeatStat_ws(x,y,te,SVMspace);
    
    % apply feature selection
    switch SVMspace
        case {'cart_HT' 'cartNoAmp_HT' 'cartNoDelay_HT'...
                'cartReal_T'...
                'polMag_T'}
            featStat = abs(featStat);
            x = x(:,featStat>prctile(featStat,10));
        case {'cart' 'cartNoAmp' 'cartNoDelay'...
                'cart_HTbSess' 'cartNoAmp_HTbSess' 'cartNoDelay_HTbSess'...
                'polMag'}
        otherwise
            error('X')
    end

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

function x = polarSpaceNormalization(x,te,SVMspace)
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

switch SVMspace
    case {'hr' 'hrNoAmp'}
        % does not apply
    case {'cart' 'cart_HT' 'cart_HTbSess'...
            'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'...
            'cartNoDelay' 'cartNoDelay_HT' 'cartNoDelay_HTbSess'...
            'cartReal_T'...
            'polMag' 'polMag_T'}
        % set to mean rho=1 and mean theta=0 in each voxel)
        switch SVMspace
            case {'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'}
                % but set rho=1 for each vector (omit any amplitude information)
                rho = 1;
            case {'cart' 'cart_HT' 'cart_HTbSess'...
                    'cartNoDelay' 'cartNoDelay_HT' 'cartNoDelay_HTbSess'...
                    'cartReal_T'...
                    'polMag' 'polMag_T'}
                rho = abs(x)./abs(mean(x(~te,:),1));
            otherwise
                error('X')
        end
        switch SVMspace
            case {'cart' 'cart_HT' 'cart_HTbSess'...
                    'cartNoAmp' 'cartNoAmp_HT' 'cartNoAmp_HTbSess'...
                    'cartReal_T'}
                theta = angle(x) - angle(mean(x(~te,:),1)); theta = wrapToPi(theta);
            case {'cartNoDelay' 'cartNoDelay_HT' 'cartNoDelay_HTbSess'...
                    'polMag' 'polMag_T'}
                % but set theta=0 for each vector (omit any delay information)
                theta = 0;
            otherwise
                error('X')
        end
        [u,v] = pol2cart(theta,rho);
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
te = false(size(x,1),1);
x = polarSpaceNormalization(x,te,SVMspace);

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