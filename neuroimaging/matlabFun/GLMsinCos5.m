function [results,dataDetrend] = GLMsinCos5(design,data,stimdur,tr,hrfmodel,hrfknobs,opt,splitIn,exclusion)

%WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING
%stimdur is internally fixed to 6 in GLMestimatemodel>fitmodel_helper at 884
%WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING

numruns = length(data);
sessLabel = cell2mat(opt.sessionLabel)';
condLabel = cat(1,ones(numruns/3,1)*1,ones(numruns/3,1)*2,ones(numruns/3,1)*3);
runOrder = cell2mat(opt.runLabel)';

sessList = unique(sessLabel);
condList = unique(condLabel);

runLabel = nan(size(sessLabel));
for sessInd = 1:length(sessList)
    for condInd = 1:length(condList)
        curInd = sessLabel==sessList(sessInd) & condLabel==condList(condInd);
        runLabel(curInd) = 1:length(runLabel(curInd));
    end
end
runList = unique(runLabel);



%% Exclude
if exist('exclusion','var')
    exclInd = exclusion.sess==sessLabel & exclusion.run==runLabel;
    
    design = design(~exclInd);
    data = data(~exclInd);
    opt.sessionLabel = opt.sessionLabel(~exclInd);
    opt.runLabel = opt.runLabel(~exclInd);
    
    
    
    
    numruns = length(data);
    sessLabel = cell2mat(opt.sessionLabel)';
    condLabel = cat(1,ones(numruns/3,1)*1,ones(numruns/3,1)*2,ones(numruns/3,1)*3);
    runOrder = cell2mat(opt.runLabel)';
    
    sessList = unique(sessLabel);
    condList = unique(condLabel);
    
    runLabel = nan(size(sessLabel));
    for sessInd = 1:length(sessList)
        for condInd = 1:length(condList)
            curInd = sessLabel==sessList(sessInd) & condLabel==condList(condInd);
            runLabel(curInd) = 1:length(runLabel(curInd));
        end
    end
    runList = unique(runLabel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEAL WITH INPUTS, ETC.

% input
if strcmp(hrfmodel,'SinCos')
    for i = 1:length(design)
        design{i} = repmat(design{i},1,2);
    end
end
if ~exist('hrfmodel','var') || isempty(hrfmodel)
    hrfmodel = 'optimize';
end
if ~exist('hrfknobs','var') || isempty(hrfknobs)
    if isequal(hrfmodel,'fir')
        hrfknobs = 20;
    elseif isequal(hrfmodel,'SinCos')
        hrfknobs = sin(2*pi*1/(stimdur*2)*(0:tr:stimdur*2-tr))';
        hrfknobs = normalizemax(hrfknobs);
    else
        hrfknobs = normalizemax(getcanonicalhrf(stimdur,tr)');
    end
end
if ~exist('opt','var') || isempty(opt)
    opt = struct();
end

% massage input
if ~iscell(design)
    design = {design};
end
if ~iscell(data)
    data = {data};
end
for p=1:length(data)
    if ~isa(data{p},'single')
        %     fprintf('*** GLMdenoisedata: converting data in run %d to single format (consider doing this before the function call to reduce memory usage). ***\n',p);
        %     data{p} = single(data{p});
        %         fprintf('*** GLMdenoisedata: ***NOT*** converting data in run %d to single format (consider doing this before the function call to reduce memory usage). ***\n',p);
        %     data{p} = single(data{p});
    end
end

% do some error checking
if any(flatten(~isfinite(data{1})))
    fprintf('*** GLMdenoisedata: WARNING: we checked the first run and found some non-finite values (e.g. NaN, Inf). unexpected results may occur due to non-finite values. please fix and re-run GLMdenoisedata. ***\n');
end

% calc
numruns = length(design);
dataclass = class(data{1});  % will always be 'single'
is3d = size(data{1},4) > 1;
if is3d
    dimdata = 3;
    dimtime = 4;
    xyzsize = sizefull(data{1},3);
else
    dimdata = 1;
    dimtime = 2;
    xyzsize = size(data{1},1);
end
numvoxels = prod(xyzsize);

% deal with defaults

if ~isfield(opt,'extraregressors') || isempty(opt.extraregressors)
    opt.extraregressors = cell(1,numruns);
end
if ~isfield(opt,'maxpolydeg') || isempty(opt.maxpolydeg)
    opt.maxpolydeg = zeros(1,numruns);
    for p=1:numruns
        opt.maxpolydeg(p) = round(((size(data{p},dimtime)*tr)/60)/2);
    end
end

if length(opt.maxpolydeg) == 1
    opt.maxpolydeg = repmat(opt.maxpolydeg,[1 numruns/splitIn]);
end
if ~iscell(opt.extraregressors)
    opt.extraregressors = {opt.extraregressors};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DETERMINE HRF

hrf = hrfknobs;
hrffitvoxels = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETRIC FITS AND ERROR ESTIMATES

for p = 1:length(design)
    ntime = size(design{p},1);                    % number of time points
    hrfknobs(:,2) = normalizemax(cos(2*pi*1/(6*2)*(0:tr:6*2-tr))');
    for i = 1:2
        tmp(:,i) = conv2(full(design{p}(:,i)),hrfknobs(:,i));  % convolve
    end
    tmp = tmp(1:ntime,:);             % extract desired subset
    cache.rawdesign{p} = tmp; clear tmp
end

if ~isfield(opt,'splitedIn')
    opt.splitedIn = 1;
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIXED-EFFECT FIT
tic
fprintf('*** FIXED-EFFECT FIT ***\n');
%% Construct design matrix
% Polynomials and motion
poly = cell(1,numruns);
polyNmotion = cell(1,numruns);
constant = cell(1,numruns);

numtime = size(data{p},dimtime)*splitIn;
for p=1:numruns/splitIn
    tmp = constructpolynomialmatrix(numtime,0:opt.maxpolydeg(p),0);
    poly((p-1)*splitIn+1:p*splitIn) = mat2cell(tmp,repmat(size(data{p},dimtime),1,splitIn),2)';
end

numtime = size(data{p},dimtime)*splitIn;
for p=1:numruns
    polyNmotion{p} = cat(2,poly{p}, ...
        opt.extraregressors{p}, ...
        []);
    tmpConstant = poly{p}; tmpConstant(:,2:end) = zeros(size(tmpConstant(:,2:end)));
    constant{p} = cat(2,tmpConstant, ...
        zeros(size(opt.extraregressors{p})), ...
        []);
end
%Reformat to remerge motion regressors that were split before
polyNmotion_b = reshape(polyNmotion,[opt.splitedIn length(polyNmotion)/opt.splitedIn])';
constant_b = reshape(constant,[opt.splitedIn length(constant)/opt.splitedIn])';
polyNmotion = cell(1,size(polyNmotion_b,1));
constant = cell(1,size(polyNmotion_b,1));
for ii = 1:size(polyNmotion_b,1)
    polyNmotion{ii} = catcell(1,polyNmotion_b(ii,:));
    constant{ii} = catcell(1,constant_b(ii,:));
end
clear polyNmotion_b constant_b

% Conditions
condDesign1 = cell(1,3);
polyNmotionDesign1 = cell(1,3);
constantDesign1 = cell(1,3);
condDesign2 = cell(1,3);
polyNmotionDesign2 = cell(1,3);
constantDesign2 = cell(1,3);
for cond = 1:3
    sess=1;
    condInd = condLabel==cond & sessLabel==sess;
    condInd_polyNmotion = condInd(1:opt.splitedIn:end);
    condDesign1{cond} = catcell(1,cache.rawdesign(condInd));
    
    polyNmotionDesign1{cond} = blkdiag(polyNmotion{condInd_polyNmotion});
    constantDesign1{cond} = blkdiag(constant{condInd_polyNmotion});
    
    
    sess=2;
    condInd = condLabel==cond & sessLabel==sess;
    condInd_polyNmotion = condInd(1:opt.splitedIn:end);
    condDesign2{cond} = catcell(1,cache.rawdesign(condInd));
    
    polyNmotionDesign2{cond} = blkdiag(polyNmotion{condInd_polyNmotion});
    constantDesign2{cond} = blkdiag(constant{condInd_polyNmotion});
end
condDesign = reshape(cat(1,condDesign1,condDesign2),1,6);
polyNmotionDesign = reshape(cat(1,polyNmotionDesign1,polyNmotionDesign2),1,6);
constantDesign = reshape(cat(1,constantDesign1,constantDesign2),1,6);


% Merge conditions, polynomials and motion
condDesign = blkdiag(condDesign{:});
% sessDesign = catcell(1,sessDesign); sessDesign(:,1:2) = [];
polyNmotionDesign = blkdiag(polyNmotionDesign{:});
constantDesign = blkdiag(constantDesign{:});

results.OLS.fixed.designmatrix = cat(2,condDesign,polyNmotionDesign);% time x regressors
% results.OLS.fixed.designmatrixCond3 = cat(2,condDesign(end*2/3+1:end,:),polyNmotionDesign(end*2/3+1:end,:));% time x regressors

% Extract some usefull params
results.OLS.fixed.designmatrixPieces.cond = cat(2,condDesign,zeros(size(polyNmotionDesign)));
results.OLS.fixed.designmatrixPieces.motion = cat(2,zeros(size(condDesign)),polyNmotionDesign);
results.OLS.fixed.designmatrixPieces.constant = cat(2,zeros(size(condDesign)),constantDesign);
%     close all
%     figure('WindowStyle','docked'); colormap gray
%     imagesc(results.OLS.fixed.designmatrix)
%     figure('WindowStyle','docked'); colormap gray
%     imagesc(results.OLS.fixed.designmatrixPieces.cond)
%     figure('WindowStyle','docked'); colormap gray
%     imagesc(results.OLS.fixed.designmatrixPieces.motion)
%     figure('WindowStyle','docked'); colormap gray
%     imagesc(results.OLS.fixed.designmatrixPieces.constant,[-1 1])


%% Estimate parameters with OLS
display('computing OLS')
results.OLS.fixed.parameters = ...
    mtimescell(olsmatrix2(results.OLS.fixed.designmatrix), ...
    cellfun(@(x) squish(x,dimdata)',data,'UniformOutput',0));  % regressors x voxels





%% Compute F session-wise
tmp = repmat(sessLabel',numtime,1);
results.OLS.fixed.designmatrixSessInd = reshape(tmp,numel(tmp),1);
tmp = repmat(condLabel',numtime,1);
results.OLS.fixed.designmatrixCondInd = reshape(tmp,numel(tmp),1);
tmp = repmat(runOrder',numtime,1);
results.OLS.fixed.designmatrixRunInd = reshape(tmp,numel(tmp),1);

for sess = 1:2
    display(['Computing F for sess' num2str(sess)])
    % Set-up full and reduced design matrices
    designMatrixFull = results.OLS.fixed.designmatrix(results.OLS.fixed.designmatrixSessInd==sess,:);
    paramInd = any(designMatrixFull,1);    
    param = results.OLS.fixed.parameters(paramInd,:);
    designMatrixFull = designMatrixFull(:,paramInd);
    designMatrixReduced = designMatrixFull;
    designMatrixReduced(:,1:6) = zeros(size(designMatrixReduced(:,1:6)));
    
    results.OLS.fixed.(['sess' num2str(sess)]).designMatrixFull = designMatrixFull;
    results.OLS.fixed.(['sess' num2str(sess)]).designMatrixReduced = designMatrixReduced;
    results.OLS.fixed.(['sess' num2str(sess)]).designMatrixConstant = results.OLS.fixed.designmatrixPieces.constant(results.OLS.fixed.designmatrixSessInd==sess,paramInd);

    
    % Compute predictions
    fullModel = zeros(numvoxels,size(designMatrixFull,1));
    reducedModel = zeros(numvoxels,size(designMatrixReduced,1));
    for vox = 1:numvoxels
        tmpModel = bsxfun(@times,param(:,vox)',designMatrixFull);
        fullModel(vox,:) = sum(tmpModel,2);
        reducedModel(vox,:) = sum(tmpModel(:,any(designMatrixReduced,1)),2);
    end
    fullModel = reshape(fullModel,  [xyzsize size(fullModel,2)]);
    reducedModel = reshape(reducedModel,  [xyzsize size(reducedModel,2)]);
   
    % Compute RSSs
    RSSfull = bsxfun(@minus, fullModel, catcell(4,data(sessLabel==sess)));
    RSSfull = bsxfun(@times, RSSfull, RSSfull);
    RSSfull = sum(RSSfull,4);
    
    RSSreduced = bsxfun(@minus, reducedModel, catcell(4,data(sessLabel==sess)));
    RSSreduced = bsxfun(@times, RSSreduced, RSSreduced);
    RSSreduced = sum(RSSreduced,4);
    
    % Compute and output F (considering all data)
    %Define degrees of freedom
    pFull = length(find(any(designMatrixFull,1)));
    pReduced = length(find(any(designMatrixReduced,1)));
    n = size(designMatrixFull,1);
    %Compute F
    results.OLS.fixed.(['sess' num2str(sess)]).F.val.F = ((RSSreduced-RSSfull)./(pFull-pReduced)) ./ (RSSfull./(n-pFull));
    results.OLS.fixed.(['sess' num2str(sess)]).F.val.F(results.OLS.fixed.(['sess' num2str(sess)]).F.val.F<=0) = eps(class(results.OLS.fixed.(['sess' num2str(sess)]).F.val.F)); % replace with smallest possible number if impossible value
    results.OLS.fixed.(['sess' num2str(sess)]).F.df.pFull = pFull;
    results.OLS.fixed.(['sess' num2str(sess)]).F.df.pReduced = pReduced;
    results.OLS.fixed.(['sess' num2str(sess)]).F.df.n = n;
    clear temp temp2 sumsq good X fullModel reducedModel RSSfull RSSreduced
    
    %% Reformat outputs
    results.OLS.fixed.(['sess' num2str(sess)]).parameters =   reshape(param',  [xyzsize size(param,1)]);
    
%     %% Convert to percent BOLD
%     results.OLS.fixed.(['sess' num2str(sess)]).constant.brain = mean(results.OLS.fixed.(['sess' num2str(sess)]).parameters(:,:,:,any(results.OLS.fixed.(['sess' num2str(sess)]).designMatrixConstant,1)).*(vectorlength(ones(numtime,1))/numtime),4); %vectorlength(ones(numtime,1))/numtime is the value in the design matrix
    results.OLS.fixed.(['sess' num2str(sess)]).constant.brain = results.OLS.fixed.(['sess' num2str(sess)]).parameters(:,:,:,any(results.OLS.fixed.(['sess' num2str(sess)]).designMatrixConstant,1)); %.*(vectorlength(ones(numtime,1))/numtime),4); %vectorlength(ones(numtime,1))/numtime is the value in the design matrix
%     con = 1./abs(results.OLS.fixed.(['sess' num2str(sess)]).constant.brain) * 100;
%     results.OLS.fixed.(['sess' num2str(sess)]).parameters = bsxfun(@times,results.OLS.fixed.(['sess' num2str(sess)]).parameters,con);
    
    %% Compute amp and delay
    [results.OLS.fixed.(['sess' num2str(sess)]).delay, results.OLS.fixed.(['sess' num2str(sess)]).amp] = cart2pol(results.OLS.fixed.(['sess' num2str(sess)]).parameters(:,:,:,1:2:6),results.OLS.fixed.(['sess' num2str(sess)]).parameters(:,:,:,2:2:6));
    toc
end

%% Detrend data
dataDetrend = [];
% display('detrending data')
% params =   reshape(results.OLS.fixed.parameters',  [xyzsize size(results.OLS.fixed.parameters,1)]);
% designmatrix = results.OLS.fixed.designmatrix;
% designmatrix(:,any(results.OLS.fixed.designmatrixPieces.cond,1)) = 0;
% confoundInd = find(any(results.OLS.fixed.designmatrixPieces.constant,1));
% dataDetrend = data;
% 
% % Detrend data
% for runInd = 1:size(data,2)
%     curRunInd = logical(results.OLS.fixed.designmatrixPieces.constant(:,confoundInd(runInd)));
%     curModel = designmatrix(curRunInd,:);
%     curData = data{runInd};
%     curDataDetrend = curData;
%     
%     for x = 1:size(curData,1)
%         for y = 1:size(curData,2)
%             for z = 1:size(curData,3)
%                 curParam = squeeze(params(x,y,z,:));
%                 constant = curParam(find(any(curModel,1),1)).*(vectorlength(ones(numtime,1))/numtime); % for %BOLD   (vectorlength(ones(numtime,1))/numtime is the value in the design matrix)
%                 curParam = repmat(curParam',size(curData,4),1);
%                 curDataDetrend(x,y,z,:) = squeeze(curData(x,y,z,:)) - sum(curParam.*curModel,2);
%                 % convert to %BOLD
%                 curDataDetrend(x,y,z,:) = curDataDetrend(x,y,z,:)./constant;
%             end
%         end
%     end
%     dataDetrend{runInd} = curDataDetrend;
% end




if isfield(opt,'sessModel')
    opt.sessModel = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE ADDITIONAL OUTPUTS

% return all the inputs (except for the data) in the output.
% also, include a new field 'datasize'.
results.inputs.design = design;
results.inputs.datasize = cellfun(@(x) size(x),data,'UniformOutput',0);
results.inputs.stimdur = stimdur;
results.inputs.tr = tr;
results.inputs.hrfmodel = hrfmodel;
results.inputs.hrfknobs = hrfknobs;
results.inputs.opt = opt;



