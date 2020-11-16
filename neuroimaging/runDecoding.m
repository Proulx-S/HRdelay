function res = runDecoding(threshType,SVMspace)
if ~exist('SVMspace','var') || isempty(SVMspace)
    SVMspace = 'cartReal'; % 'cart', 'cartReal', 'cartRealFixedDelay', 'cartImag', 'pol', 'polMag' or 'polDelay'
end
if ~exist('threshType','var') || isempty(threshType)
    threshType = 'fdr'; % 'none', 'p' or 'fdr'
end
threshVal = 0.05;


repoPath = 'C:\Users\sebas\OneDrive - McGill University\dataBig';
dataDir = 'C-derived\DecodingHR';
funPath = fullfile(repoPath,dataDir,'fun');
funLevel_in = 'zSin';
funLevel_out = 'zSin/decoding';
subjList = {'02jp' '03sk' '04sp' '05bm' '06sb' '07bj'}';
fileSuffix_in = '_maskSinAndHrFit.mat';
fileSuffix_out = '_decoding.mat';



disp(['IN: Sinusoidal BOLD responses from anatomical V1 ROI (' fullfile(dataDir,funLevel_in) ')'])
disp('F(IN)=OUT: threshold included voxels and decode ROI response pattern to predict stimulus orientation')
disp(['OUT: figures and stats (' fullfile(dataDir,funLevel_in) ')'])



%% Load data
dAll = cell(size(subjList,1),1);
for subjInd = 1:size(subjList,1)
    load(fullfile(funPath,funLevel_in,[subjList{subjInd} fileSuffix_in]),'d');
    dAll{subjInd} = d;
end
d = dAll; clear dAll

%% Pipe data
% Threshold and average voxels in cartesian space
dP = d;
dP2 = cell(size(dP));
for subjInd = 1:length(dP)
    switch threshType
        % no voxel selection
        case 'none'
            % voxel selection cross-validated between sessions
            indSess1 = true(size(dP{subjInd}.sess1.F));
            indSess2 = true(size(dP{subjInd}.sess1.F));
        case 'p'
            indSess1 = dP{subjInd}.sess1.P<threshVal;
            indSess2 = dP{subjInd}.sess2.P<threshVal;
        case 'fdr'
            indSess1 = dP{subjInd}.sess1.FDR<threshVal;
            indSess2 = dP{subjInd}.sess2.FDR<threshVal;
        otherwise
            error('X')
    end
    dP{subjInd}.sess1.xData = dP{subjInd}.sess1.xData(:,indSess2,:);
    dP{subjInd}.sess1.F     = dP{subjInd}.sess1.F(:,indSess2,:);
    dP{subjInd}.sess1.FDR   = dP{subjInd}.sess1.FDR(:,indSess2,:);
    dP{subjInd}.sess1.P     = dP{subjInd}.sess1.P(:,indSess2,:);
    dP{subjInd}.sess2.xData = dP{subjInd}.sess2.xData(:,indSess1,:);
    dP{subjInd}.sess2.F     = dP{subjInd}.sess2.F(:,indSess1,:);
    dP{subjInd}.sess2.FDR   = dP{subjInd}.sess2.FDR(:,indSess1,:);
    dP{subjInd}.sess2.P     = dP{subjInd}.sess2.P(:,indSess1,:);
    
    sessList = fields(dP{subjInd});
    for sessInd = 1:length(sessList)
        dP2{subjInd,sessInd} = dP{subjInd}.(sessList{sessInd});
    end
    dP{subjInd} = [];
end
dP = dP2; clear dP2



%% Run svm
res.nVox = nan(size(dP));
res.nDim = nan(size(dP));
res.acc  = nan(size(dP));
res.nObs = nan(size(dP));
res.p = nan(size(dP));
res.subjList = subjList;
for i = 1:numel(dP)
    nSamplePaired = size(dP{i}.xData,1);
    
    x1 = dP{i}.xData(:,:,1);
    y1 = 1.*ones(nSamplePaired,1);
    k1 = (1:nSamplePaired)';
    
    x2 = dP{i}.xData(:,:,2);
    y2 = 2.*ones(nSamplePaired,1);
    k2 = (1:nSamplePaired)';
    
    
    y = cat(1,y1,y2); clear y1 y2
    k = cat(1,k1,k2); clear k1 k2
    
%     %add an artificial offset for validation purposes
%     x2std = mean(abs(x2(:)))./std(abs(x2(:)));
%     [X,Y] = pol2cart(angle(x2),abs(x2) + x2std*0.2);
%     x2 = complex(X,Y);
      

    % SVM
    kList = unique(k);
    yTr = nan(length(y),length(kList));
    yTe = nan(length(y),1);
    yHatTr = nan(length(y),length(kList));
    yHatTe = nan(length(y),1);
    for kInd = 1:length(kList)
        x = cat(1,x1,x2);
        
        % split train and test
        te = k==kList(kInd);
        
        % polar space normalization (rho=1, theta=0)
        if ~strcmp(SVMspace,'cartRealFixedDelay')
            rho = abs(mean(x(~te,:),1));
            theta = angle(mean(x(~te,:),1));
        else
            rho = abs(mean(mean(x(~te,:),1),2));
            theta = angle(mean(mean(x(~te,:),1),2));
        end
        [X,Y] = pol2cart(angle(x)-theta,abs(x)./rho);
        x = complex(X,Y); clear X Y
        
        switch SVMspace
            case 'cart'
                x = x./std(x(~te,:),[],1) - mean(x(~te,:),1);
                x = cat(2,real(x),imag(x));
            case {'cartReal'  'cartRealFixedDelay'}
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
        model = svmtrain(y(~te,:),x(~te,:),'-t 2 -q');
%         w = model.sv_coef'*model.SVs;
%         b = model.rho;
%         yHat = cat(2,real(x(~te,:)),imag(x(~te,:)))*w';
        [yTr(~te,kInd), ~, yHatTr(~te,kInd)] = svmpredict(y(~te,:),x(~te,:),model,'-q');
        [yTe(te,1), ~, yHatTe(te,1)] = svmpredict(y(te,:),x(te,:),model,'-q');
    end
    res.nVox(i) = size(x1,2);
    res.nDim(i) = size(x,2);
    res.nObs(i) = length(y);
    res.acc(i) = sum(yTe==y)./res.nObs(i);
    res.p(i) = binocdf(res.acc(i).*res.nObs(i),res.nObs(i),0.5,'upper');
end
hit = sum(res.acc(:).*res.nObs(:));
n = sum(res.nObs(:));
p = binocdf(sum(res.acc(:).*res.nObs(:)),sum(res.nObs(:)),0.5,'upper');
disp('---');
disp(['SVM space: ' SVMspace '; Vox selection: ' threshType]);
disp(['Group accuracy = ' num2str(hit) '/' num2str(n) ' (' num2str(hit/n*100,'%0.1f') '%; binomial p=' num2str(p,'%0.3f') ')'])
