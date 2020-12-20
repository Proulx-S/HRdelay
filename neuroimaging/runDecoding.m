function res = runDecoding(SVMspace,nPerm,figOption)
if ~exist('figOption','var') || isempty(figOption)
    figOption.save = 0;
    figOption.subj = 1; % 'all' or subjInd
end
if ~exist('SVMspace','var') || isempty(SVMspace)
    SVMspace = 'cart_HT'; % 'hr' 'hrNoAmp' 'cart' 'cartNoAmp' cartNoAmp_HT 'cartReal', 'cartImag', 'pol', 'polMag' 'polMag_T' or 'polDelay'
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
fileSuffix = '_defineAndShowMasks.mat';

%make sure everything is forward slash for mac, linux pc compatibility
for tmpPath = {'repoPath' 'dataDir' 'funPath'}
    eval([char(tmpPath) '(strfind(' char(tmpPath) ',''\''))=''/'';']);
end


%% Preload param
tmp = dir(fullfile(funPath,funLevel,'*.mat'));
for i = 1:length(tmp)
    load(fullfile(funPath,funLevel,tmp(i).name),'param')
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
    load(fullfile(funPath,funLevel,[subjList{subjInd} fileSuffix]),'dC');
    dCAll{subjInd} = dC;
end
dC = dCAll; clear dAll

% save('C:\Users\sebas\Desktop\main','dC')


%% Further between-session feature selection
% all managed in previous steps
dP = cell(length(dC),2);
for subjInd = 1:length(dC)
    switch SVMspace
        case {'cart' 'cart_HT' 'cartNoAmp' 'cartNoAmp_HT' 'polMag' 'polMag_T'}
            % Does not apply
            %Between-session feature selection of voxels activated in any
            %condition 'anyCondActivation' is currently done in previous steps:
%             featSelContrast1 = param.featSelContrast1.name;
%             threshType = param.featSelContrast1.threshType;
%             threshVal = param.featSelContrast1.threshVal;
            %It should however be done here to allow other kinds of
            %between-subject feature selection (e.g. 'cart_HTbSess' below)
            sessList = fields(dC{subjInd});
            for sessInd = 1:length(sessList)
                dP{subjInd,sessInd} = dC{subjInd}.(sessList{sessInd});
            end
            dC{subjInd} = [];
            
        case 'cart_HTbSess'
            error('code not complete')
            x = dC{subjInd}.(sessList{sessInd}).data;
            y = []; for condInd = 1:size(x,3); y = cat(3,y,ones(size(x,1),1).*condInd); end
            x2 = []; y2 = [];
            for condInd = 1:size(x,3)
                x2 = cat(1,x2,x(:,:,condInd));
                y2 = cat(1,y2,y(:,:,condInd));
            end
            x = x2; y = y2; clear x2 y2
            featStat = nan(1,size(x,2));
            for voxInd = 1:size(x,2)
                xTmp = cat(1,x(y==1,voxInd),x(y==2,voxInd));
                stats = T2Hot2d([real(xTmp) imag(xTmp)]);
                featStat(voxInd) = stats.T2;
            end
        otherwise
            error('X')
    end
end
clear dC

%% Example plot of trigonometric (polar representation)
switch param.featSelContrast1.name
    case 'anyCondActivation'
        [~,b] = max(dP{1}.([param.featSelContrast1.name '__F']));
    otherwise
        error('X')
end
x = dP{1}.data(:,b,1:3);

switch SVMspace
    case {'hr' 'hrNoAmp'}
        % does not apply
    case {'cart' 'cart_HT' 'cartNoAmp' 'cartNoAmp_HT' 'polMag' 'polMag_T'}
        % set to mean rho=1 and mean theta=0 in each voxel)
        switch SVMspace
            case {'cartNoAmp' 'cartNoAmp_HT'}
                % but set rho=1 for each vector (omit any amplitude information)
                rho = 1;
            case {'cart' 'cart_HT' 'polMag' 'polMag_T'}
                rho = abs(x)./abs(mean(x,1));
            otherwise
                error('X')
        end
        switch SVMspace
            case {'cart' 'cart_HT' 'cartNoAmp' 'cartNoAmp_HT'}
                % but set rho=1 for each vector (omit any amplitude information)
                theta = angle(x) - angle(mean(x,1)); theta = wrapToPi(theta);
            case {'polMag' 'polMag_T'}
                theta = 0;
            otherwise
                error('X')
        end
        [u,v] = pol2cart(theta,rho);
        xAfter = complex(u,v); clear u v
    otherwise
        error('X')
end

f = figure('WindowStyle','docked');
subplot(1,2,1); clear hPP
% polarplot(angle(x(:)),abs(x(:)),'.'); hold on
for condInd = 1:3
    hPP(condInd) = polarplot(angle(x(:,:,condInd)),abs(x(:,:,condInd)),'o'); hold on
    hPP(condInd).MarkerFaceColor = hPP(condInd).Color;
    hPP(condInd).MarkerEdgeColor = 'w';
    hPP(condInd).MarkerSize = 3.5;
end

ax1 = gca;
subplot(1,2,2);
% polarplot(angle(xAfter(:)),abs(xAfter(:)),'.'); hold on
for condInd = 1:3
    hPP(condInd) = polarplot(angle(xAfter(:,:,condInd)),abs(xAfter(:,:,condInd)),'o'); hold on
    hPP(condInd).MarkerFaceColor = hPP(condInd).Color;
    hPP(condInd).MarkerEdgeColor = 'w';
    hPP(condInd).MarkerSize = 3.5;
end
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
ax.ThetaTickLabel = 12-ax.ThetaTick(1:end)/360*12;
ax.ThetaTickLabel(1,:) = '0 ';
% ax.ThetaAxis.Label.String = {'delay' '(sec)'};
% ax.ThetaAxis.Label.Rotation = 0;
% ax.ThetaAxis.Label.HorizontalAlignment = 'left';
ax.RAxis.Label.String = 'amp (%BOLD)';
ax.RAxis.Label.Rotation = 80;
ax.Title.String = 'after';

hSup = suptitle({'Polar space normalization' SVMspace});
hSup.Interpreter = 'none';
drawnow

if figOption.save
    filename = fullfile(pwd,mfilename);
    if ~exist(filename,'dir'); mkdir(filename); end
    filename = fullfile(filename,SVMspace);
    f.Color = 'none';
    set(findobj(f.Children,'type','Axes'),'color','none')
    saveas(f,[filename '.svg']); disp([filename '.svg'])
    f.Color = 'w';
    set(findobj(f.Children,'type','Axes'),'color','w')
    saveas(f,filename); disp([filename '.fig'])
    saveas(f,filename); disp([filename '.jpg'])
end

%% Some more independant-sample normalization
switch SVMspace
    case 'hrNoAmp' % scale each hr to 1 accroding to sin fit
        for i = 1:numel(dP)
            dP{i}.hr = dP{i}.hr./abs(dP{i}.data./100);
        end
    case {'cart' 'cart_HT' 'cartNoAmp' 'cartNoAmp_HT' 'polMag' 'polMag_T'}
        % Does not apply
    otherwise
end


%% Run svm
if ~doPerm
    res.nVox = nan(size(dP));
    res.nDim = nan(size(dP));
    res.acc  = nan(size(dP));
    res.nObs = nan(size(dP));
    res.p = nan(size(dP));
    res.subjList = subjList;
else
    res.perm.acc = nan([nPerm size(res.acc)]);
end
for i = 1:numel(dP)
    if doPerm
        disp(['for sess ' num2str(i) ' of ' num2str(numel(dP))])
        tic
    end
    % Define x(data), y(label) and k(xValFolds)
    switch SVMspace
        case {'hr' 'hrNoAmp'}
            nSamplePaired = size(dP{i}.hr,1);
            x1 = dP{i}.hr(:,:,1,:); x1 = x1(:,:);
            x2 = dP{i}.hr(:,:,2,:); x2 = x2(:,:);
        case {'cart' 'cart_HT' 'cartNoAmp' 'cartNoAmp_HT' 'polMag' 'polMag_T'}
            nSamplePaired = size(dP{i}.data,1);
            x1 = dP{i}.data(:,:,1);
            x2 = dP{i}.data(:,:,2);
        otherwise
            error('X')
    end
    y1 = 1.*ones(nSamplePaired,1);
    k1 = (1:nSamplePaired)';
    y2 = 2.*ones(nSamplePaired,1);
    k2 = (1:nSamplePaired)';
    
    x = cat(1,x1,x2); clear x1 x2
    y = cat(1,y1,y2); clear y1 y2
    k = cat(1,k1,k2); clear k1 k2
    
    % SVM
    if ~doPerm
        %cross-validated SVM
        [yTe,d] = xValSVM(x,y,k,SVMspace);
    else
        % with permutations
        yTe = nan(length(y),nPerm);
        y1 = y(y==1);
        y2 = y(y==2);
        parfor permInd = 1:nPerm
            %permute labels
            y = cat(2,y1,y2);
            for sInd = 1:size(y,1); y(sInd,:) = y(sInd,randperm(2)); end
            y = cat(1,y(:,1),y(:,2));
            %cross-validated SVM
            yTe(:,permInd) = xValSVM(x,y,k,SVMspace);
        end
        y = cat(1,y1,y2); clear y1 y2
    end
    
    % Evaluate SVM
    if ~doPerm
        res.nObs(i) = length(y);
        res.acc(i) = sum(yTe==y)./res.nObs(i);
        res.p(i) = binocdf(res.acc(i).*res.nObs(i),res.nObs(i),0.5,'upper');
        res.nDim(i) = mean(d);
        switch SVMspace
            case 'hr'
                res.nVox1(i) = size(dP{i}.hr,2);
                res.nVox2(i) = mean(d)./size(dP{i}.hr,4);
            case {'cart' 'cart_HT' 'cartNoAmp' 'cartNoAmp_HT'}
                res.nVox1(i) = size(dP{i}.data,2);
                res.nVox2(i) = mean(d)./2;
            case {'polMag' 'polMag_T'}
                res.nVox1(i) = size(dP{i}.data,2);
                res.nVox2(i) = mean(d);
            otherwise
                error('X')
        end
    else
        res.perm.acc(:,i) = sum(yTe==y,1)./res.nObs(i);
    end
    if doPerm
        toc
    end
end
%% Add info
if ~doPerm
    res.info = 'subj x sess';
    res.summary.SVMspace = SVMspace;
    res.summary.hit = sum(res.acc(:).*res.nObs(:));
    res.summary.nObs = sum(res.nObs(:));
    res.summary.acc = res.summary.hit/res.summary.nObs;
    [~,pci] = binofit(res.summary.nObs/2,res.summary.nObs,0.1);
    res.summary.accThresh = pci(2);
    res.summary.p = binocdf(res.summary.hit,res.summary.nObs,0.5,'upper');
    disp('Group results:')
    disp(['  hit    =' num2str(res.summary.hit) '/' num2str(res.summary.nObs)])
    disp(['  acc    =' num2str(res.summary.acc*100,'%0.2f%%')])
    disp(' binomial stats')
    disp(['  thresh =' num2str(res.summary.accThresh*100,'%0.2f%%')])
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


function [yTe,d] = xValSVM(x,y,k,SVMspace)
X = x;
kList = unique(k);
% yTr = nan(length(y),length(kList));
yTe = nan(length(y),1);
% yHatTr = nan(length(y),length(kList));
yHatTe = nan(length(y),1);
d = nan(length(kList),1);
for kInd = 1:length(kList)
    x = X;
    
    % Split train and test
    te = k==kList(kInd);
    
    % Polar space normalization
    switch SVMspace
        case {'hr' 'hrNoAmp'}
            % does not apply
        case {'cart' 'cart_HT' 'cartNoAmp' 'cartNoAmp_HT' 'polMag' 'polMag_T'}
            % set to mean rho=1 and mean theta=0 in each voxel (normalize)
            switch SVMspace
                case {'cartNoAmp' 'cartNoAmp_HT'}
                    % but set rho=1 for each vector (omit any amplitude information)
                    rho = 1;
                case {'cart' 'cart_HT'...
                        'cartNoDelay' 'cartNoDelay_HT'...
                        'polMag' 'polMag_T'}
                    rho = abs(x)./abs(mean(x(~te,:),1));
                otherwise
                    error('X')
            end
            switch SVMspace
                case {'cartNoDelay' 'cartNoDelay_HT'...
                        'polMag' 'polMag_T'}
                    % but set theta=0 for each vector (omit any delay information)
                    theta = 0;
                case {'cart' 'cart_HT'...
                        'cartNoAmp' 'cartNoAmp_HT'}
                    theta = angle(x) - angle(mean(x(~te,:),1)); theta = wrapToPi(theta);
                otherwise
                    error('X')
            end
            [u,v] = pol2cart(theta,rho);
            x = complex(u,v); clear u v
        otherwise
            error('X')
    end
    
    % Within-session feature selection
    % get feature selection stats
    switch SVMspace
        case 'cart_HT'
            featStat = nan(1,size(x,2));
            for voxInd = 1:size(x,2)
                xTmp = cat(1,x(~te & y==1,voxInd),x(~te & y==2,voxInd));
%                 stats = T2Hot2d([zscore(real(xTmp)) zscore(imag(xTmp))]);
                stats = T2Hot2d([real(xTmp) imag(xTmp)]);
                featStat(voxInd) = stats.T2;
            end            
        case 'cartNoAmp_HT'
            featStat = nan(1,size(x,2));
            for voxInd = 1:size(x,2)
                [~, F] = circ_htest(angle(x(~te & y==1,voxInd)), angle(x(~te & y==2,voxInd)));
                featStat(voxInd) = F;
            end
        case {'cartNoDelay_HT' 'polMag_T'}
            [~,~,~,STATS] = ttest(abs(x(~te & y==1,:)),abs(x(~te & y==2,:)));
            featStat = abs(STATS.tstat);
        case {'cart' 'cartNoAmp' 'cartNoDelay' 'polMag'}
        otherwise
            error('X')
    end
    % apply feature selection
    switch SVMspace
        case {'cart_HT' 'cartNoAmp_HT' 'polMag_T'} % further feature selection
            x = x(:,featStat>prctile(featStat,10));
        case {'cart' 'cartNoAmp' 'polMag'}
        otherwise
            error('X')
    end
    
    % Cocktail bank normalization
    switch SVMspace
        case {'hr' 'hrNoAmp'}
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
        case {'cart' 'cart_HT' 'cartNoAmp' 'cartNoAmp_HT'}
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
            x = cat(2,real(x),imag(x));
        case 'cartReal'
            x = real(x);
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
        case 'cartImag'
            x = imag(x);
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
        case 'pol'
            x = cat(2,angle(x),abs(x));
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
        case {'polMag' 'polMag_T'}
            x = abs(x);
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
        case 'polDelay'
            x = angle(x);
            x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
        otherwise
            error('x')
    end
    
    % RunSVM
    model = svmtrain(y(~te,:),x(~te,:),'-t 0 -q');
%     w = model.sv_coef'*model.SVs;
%     b = model.rho;
%     yHat = cat(2,real(x(~te,:)),imag(x(~te,:)))*w';

%     [yTr(~te,kInd), ~, yHatTr(~te,kInd)] = svmpredict(y(~te,:),x(~te,:),model,'-q');
    [yTe(te,1), ~, yHatTe(te,1)] = svmpredict(y(te,:),x(te,:),model,'-q');
    d(kInd) = size(x,2);
end