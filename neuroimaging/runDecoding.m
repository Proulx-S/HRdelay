function res = runDecoding(SVMspace)
close all
if ~exist('SVMspace','var') || isempty(SVMspace)
    SVMspace = 'cartNoAmp_HT'; % 'hr' 'hrNoAmp' 'cart' 'cartNoAmp' cartNoAmp_HT 'cartReal', 'cartImag', 'pol', 'polMag' or 'polDelay'
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


%% Further feature selection
%not now
dP = cell(length(dC),2);
for subjInd = 1:length(dC)
    sessList = fields(dC{subjInd});
    for sessInd = 1:length(sessList)
        dP{subjInd,sessInd} = dC{subjInd}.(sessList{sessInd});
    end
    dC{subjInd} = [];
end
clear dC
% %% Pipe data
% % Threshold and average voxels in cartesian space
% dP = d;
% dP2 = cell(size(dP));
% for subjInd = 1:length(dP)
%     for sessInd = 1:2
%         switch featSelType
%             % voxel selection cross-validated between sessions
%             case 'none' % no voxel selection
%                 indTmp = true(size(dP{subjInd}.(['sess' num2str(sessInd)]).F));
%             case 'respF_p'
%                 indTmp = dP{subjInd}.(['sess' num2str(sessInd)]).P<threshVal;
%             case 'respF_fdr'
%                 indTmp = dP{subjInd}.(['sess' num2str(sessInd)]).FDR<threshVal;
%             case 'oriT_p'
%                 error('codeThat')
%             case 'T2_nVoxAsF'
%                 switch SVMspace
%                     case {'cart' 'cartReal' 'cartImag'}
%                         sz = size(dP{subjInd}.(['sess' num2str(sessInd)]).xData(:,:,1:2));
%                         stats = nan(1,sz(2));
%                         %                         y = zeros([sz(1) 1 sz(3)]);
%                         %                         y(:,:,1) = 1; y(:,:,2) = 2;
%                         parfor voxInd = 1:sz(2)
%                             x = cat(2,real(dP{subjInd}.(['sess' num2str(sessInd)]).xData(:,voxInd,1:2)),...
%                                 imag(dP{subjInd}.(['sess' num2str(sessInd)]).xData(:,voxInd,1:2)));
%                             %                         x = cat(1,cat(2,y(:,:,1),x(:,:,1)),cat(2,y(:,:,2),x(:,:,2)));
%                             %                         stats = T2Hot2iho(x);
%                             x = cat(1,x(:,:,1),x(:,:,2));
%                             tmp = T2Hot2d(x);
%                             stats(voxInd) = tmp.T2;
%                         end
%                     case {'pol' 'polMag' 'polDelay'}
%                         sz = size(dP{subjInd}.(['sess' num2str(sessInd)]).xData(:,:,1:2));
%                         stats = nan(1,sz(2));
%                         x = dP{subjInd}.(['sess' num2str(sessInd)]).xData(:,:,1:2);
%                         parfor voxInd = 1:sz(2)
%                             x1 = abs(x(:,voxInd,:));
%                             x2 = angle(x(:,voxInd,:)); x2 = wrapToPi(x2-mean(x2(:)));
%                             xT2 = cat(2,x1,x2);
%                             xT2 = cat(1,xT2(:,:,1),xT2(:,:,2));
%                             tmp = T2Hot2d(xT2);
%                             stats(voxInd) = tmp.T2;
%                         end
%                     otherwise
%                         error('x')
%                 end
%                 indTmp = false(size(dP{subjInd}.(['sess' num2str(sessInd)]).F));
%                 [~,b] = sort(abs(stats),'descend');
%                 indTmp(b(1:sum(dP{subjInd}.(['sess' num2str(sessInd)]).FDR<threshVal))) = true;
%             case 'oriT_nVoxAsF'
%                 switch SVMspace
%                     case 'cart'
%                         error('codeThat')
%                         % should use HotellingT2.m
%                     case 'cartReal'
%                         x = real(dP{subjInd}.(['sess' num2str(sessInd)]).xData(:,:,1:2));
%                         [~,~,~,STATS] = ttest(x(:,:,1),x(:,:,2));
%                         stats = STATS.tstat;
%                     case 'cartRealFixedDelay'
%                         error('codeThat')
%                     case 'cartImag'
%                         x = imag(dP{subjInd}.(['sess' num2str(sessInd)]).xData(:,:,1:2));
%                         [~,~,~,STATS] = ttest(x(:,:,1),x(:,:,2));
%                         stats = STATS.tstat;
%                     case 'pol'
%                         % should use HotellingT2.m instead
% %                         x = dP{subjInd}.(['sess' num2str(sessInd)]).xData(:,:,1:2);
% %                         statsCirc = nan([1 size(x,2)]);
% %                         for voxInd = 1:size(x,2)
% %                             [~, statsCirc(voxInd)] = circ_htest(angle(x(:,1,1)),angle(x(:,1,2)));
% %                         end
% %                         [~,~,~,STATS] = ttest(abs(x(:,:,1)),abs(x(:,:,2)));
% %                         stats = STATS.tstat;
% %                         stats = cat(2,statsCirc,stats);
%                     case 'polMag'
%                         x = abs(dP{subjInd}.(['sess' num2str(sessInd)]).xData(:,:,1:2));
%                         [~,~,~,STATS] = ttest(x(:,:,1),x(:,:,2));
%                         stats = STATS.tstat;
%                     case 'polDelay'
%                         x = angle(dP{subjInd}.(['sess' num2str(sessInd)]).xData(:,:,1:2));
%                         statsCirc = nan([1 size(x,2)]);
%                         for voxInd = 1:size(x,2)
%                             [~, statsCirc(voxInd)] = circ_htest(angle(x(:,voxInd,1)),angle(x(:,voxInd,2)));
%                         end
%                         stats = statsCirc;
%                     otherwise
%                         error('x')
%                 end
%                 indTmp = false(size(dP{subjInd}.(['sess' num2str(sessInd)]).F));
%                 [~,b] = sort(abs(stats),'descend');
%                 indTmp(b(1:sum(dP{subjInd}.(['sess' num2str(sessInd)]).FDR<threshVal))) = true;
%             otherwise
%                 error('X')
%         end
%         eval(['indSess' num2str(sessInd) ' = indTmp;']);
%     end
%     dP{subjInd}.sess1.xData = dP{subjInd}.sess1.xData(:,indSess2,:);
%     dP{subjInd}.sess1.F     = dP{subjInd}.sess1.F(:,indSess2,:);
%     dP{subjInd}.sess1.FDR   = dP{subjInd}.sess1.FDR(:,indSess2,:);
%     dP{subjInd}.sess1.P     = dP{subjInd}.sess1.P(:,indSess2,:);
%     dP{subjInd}.sess2.xData = dP{subjInd}.sess2.xData(:,indSess1,:);
%     dP{subjInd}.sess2.F     = dP{subjInd}.sess2.F(:,indSess1,:);
%     dP{subjInd}.sess2.FDR   = dP{subjInd}.sess2.FDR(:,indSess1,:);
%     dP{subjInd}.sess2.P     = dP{subjInd}.sess2.P(:,indSess1,:);
%     
%     sessList = fields(dP{subjInd});
%     for sessInd = 1:length(sessList)
%         dP2{subjInd,sessInd} = dP{subjInd}.(sessList{sessInd});
%     end
%     dP{subjInd} = [];
% end
% dP = dP2; clear dP2



% Scale each hr to 1 accroding to sin fit
switch SVMspace
    case 'hrNoAmp'
        for i = 1:numel(dP)
            dP{i}.hr = dP{i}.hr./abs(dP{i}.data./100);
        end
    case {'cart' 'cart_HT' 'cartNoAmp' 'cartNoAmp_HT'}
        % Does not apply
    otherwise
end

% % remove global delay (because trigger might have screwed up)
% switch SVMspace
%     case {'hr' 'hrNoAmp'}
%     otherwise
%         for i = 1:numel(dP)
%             theta = wrapToPi(angle(dP{i}.data) - angle(mean(dP{i}.data,2)));
%             rho = abs(dP{i}.data);
%             [X,Y] = pol2cart(theta,rho);
%             dP{i}.data = complex(X,Y); clear X Y
%         end
% end



%% Run svm
res.nVox = nan(size(dP));
res.nDim = nan(size(dP));
res.acc  = nan(size(dP));
res.nObs = nan(size(dP));
res.p = nan(size(dP));
res.subjList = subjList;
for i = 1:numel(dP)
    % Define x(data), y(label) and k(xValFolds)
    switch SVMspace
        case {'hr' 'hrNoAmp'}
            nSamplePaired = size(dP{i}.hr,1);
            x1 = dP{i}.hr(:,:,1,:); x1 = x1(:,:);
            x2 = dP{i}.hr(:,:,2,:); x2 = x2(:,:);
        case {'cart' 'cart_HT' 'cartNoAmp' 'cartNoAmp_HT'}
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
    
    y = cat(1,y1,y2); clear y1 y2
    k = cat(1,k1,k2); clear k1 k2
        
    % SVM
    kList = unique(k);
    yTr = nan(length(y),length(kList));
    yTe = nan(length(y),1);
    yHatTr = nan(length(y),length(kList));
    yHatTe = nan(length(y),1);
    for kInd = 1:length(kList)
        x = cat(1,x1,x2);
%         hPP = polarplot(x(:),'.')
        
        % split train and test
        te = k==kList(kInd);
        
        
                
        % Within-session feature selection
        switch SVMspace
            case {'cart_HT' 'cartNoAmp_HT'} % further feature selection
                featStat = nan(1,size(x,2));
                for voxInd = 1:size(x,2)
                    xTmp = cat(1,x(~te & y==1,voxInd),x(~te & y==2,voxInd));
                    stats = T2Hot2d([real(xTmp) imag(xTmp)]);
                    featStat(voxInd) = stats.T2;
                end
                [~,b] = sort(featStat,'descend');
                x = x(:,b(1:round(end*0.9)));
            otherwise
                error('X')
        end
        
        
        % Polar space normalization
        switch SVMspace
            case {'hr' 'hrNoAmp'}
                % does not apply
            case {'cart' 'cartNoAmp' 'cartNoAmp_HT'}
                % set to mean rho=1 and mean theta=0 in each voxel)
                switch SVMspace
                    case {'cartNoAmp' 'cartNoAmp_HT'} % but set rho=1 for each vector (omit any amplitude information)
                        rho = 1;
                        theta = angle(x) - angle(mean(x(~te,:),1)); theta = wrapToPi(theta);
                    case {'cart' 'cart_HT'}
                        rho = abs(x)./abs(mean(x(~te,:),1));
                        theta = angle(x) - angle(mean(x(~te,:),1)); theta = wrapToPi(theta);
                    otherwise
                        error('X')
                end
                [X,Y] = pol2cart(theta,rho);
                x = complex(X,Y); clear X Y
            otherwise
                error('X')
        end

        % cocktail bank normalization
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
            case 'polMag'
                x = abs(x);
                x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
            case 'polDelay'
                x = angle(x);
                x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
            otherwise
                error('x')
        end
        
        % runSVM
        model = svmtrain(y(~te,:),x(~te,:),'-t 0 -q');
%         w = model.sv_coef'*model.SVs;
%         b = model.rho;
%         yHat = cat(2,real(x(~te,:)),imag(x(~te,:)))*w';
        [yTr(~te,kInd), ~, yHatTr(~te,kInd)] = svmpredict(y(~te,:),x(~te,:),model,'-q');
        [yTe(te,1), ~, yHatTe(te,1)] = svmpredict(y(te,:),x(te,:),model,'-q');
    end
    switch SVMspace
        case 'hr'
            res.nVox(i) = size(x1,2)./size(dP{1,1}.hr,4);
        case {'cart' 'cart_HT' 'cartNoAmp' 'cartNoAmp_HT'}
            res.nVox(i) = size(x1,2);
        otherwise
            error('X')
    end
    res.nDim(i) = size(x,2);
    res.nObs(i) = length(y);
    res.acc(i) = sum(yTe==y)./res.nObs(i);
    res.p(i) = binocdf(res.acc(i).*res.nObs(i),res.nObs(i),0.5,'upper');
end
hit = sum(res.acc(:).*res.nObs(:));
n = sum(res.nObs(:));
p = binocdf(sum(res.acc(:).*res.nObs(:)),sum(res.nObs(:)),0.5,'upper');
disp('---');
disp(['SVM space: ' SVMspace]);
disp(['Group accuracy = ' num2str(hit) '/' num2str(n) ' (' num2str(hit/n*100,'%0.1f') '%; binomial p=' num2str(p,'%0.3f') ')'])
