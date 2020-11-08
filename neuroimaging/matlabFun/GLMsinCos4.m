function [results,detrendData] = GLMsinCos4(design,data,stimdur,tr,hrfmodel,hrfknobs,opt)

%WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING
%stimdur is internally fixed to 6 in GLMestimatemodel>fitmodel_helper at 884
%WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING


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
end;
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
    opt.maxpolydeg = repmat(opt.maxpolydeg,[1 numruns]);
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





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MIXED-EFFECT FIT
tic
fprintf('*** MIXED-EFFECT FIT ***\n');
%% Construct design matrix
polyNmotion = cell(1,numruns);
constant = cell(1,numruns);
for p=1:numruns
    numtime = size(data{p},dimtime);
    polyNmotion{p} = cat(2,constructpolynomialmatrix(numtime,0:opt.maxpolydeg(p)), ...
        opt.extraregressors{p}, ...
        []);
    tmpConstant = constructpolynomialmatrix(numtime,0:opt.maxpolydeg(p)); tmpConstant(:,2:end) = zeros(size(tmpConstant(:,2:end)));
    constant{p} = cat(2,tmpConstant, ...
        zeros(size(opt.extraregressors{p})), ...
        []);
end
% Reformat to remerge motion regressors that were split before
polyNmotion_b = reshape(polyNmotion,[opt.splitedIn length(polyNmotion)/opt.splitedIn])';
constant_b = reshape(constant,[opt.splitedIn length(constant)/opt.splitedIn])';
polyNmotion = cell(1,size(polyNmotion_b,1));
constant = cell(1,size(polyNmotion_b,1));
for ii = 1:size(polyNmotion_b,1)
    polyNmotion{ii} = catcell(1,polyNmotion_b(ii,:));
    constant{ii} = catcell(1,constant_b(ii,:));
end
clear polyNmotion_b constant_b

% Conditions and Sessions
condLabel = cat(1,ones(numruns/3,1)*1,ones(numruns/3,1)*2,ones(numruns/3,1)*3);

condDesign = cell(1,3);
polyNmotionDesign = cell(1,3);
constantDesign = cell(1,3);
for cond = 1:3
    condInd = condLabel==cond;
    condInd_polyNmotion = condInd(1:opt.splitedIn:end);
    condDesign{cond} = blkdiag(cache.rawdesign{condInd});
    
    polyNmotionDesign{cond} = blkdiag(polyNmotion{condInd_polyNmotion});
    constantDesign{cond} = blkdiag(constant{condInd_polyNmotion});
end

% Merge conditions, sessions, polynomials and motion
condDesign = blkdiag(condDesign{:});
polyNmotionDesign = blkdiag(polyNmotionDesign{:});
constantDesign = blkdiag(constantDesign{:});

results.OLS.mixed.designmatrix = cat(2,condDesign,polyNmotionDesign);% time x regressors

% Extract some usefull params
results.OLS.mixed.designmatrixPieces.cond = cat(2,condDesign,zeros(size(polyNmotionDesign)));
results.OLS.mixed.designmatrixPieces.motion = cat(2,zeros(size(condDesign)),polyNmotionDesign);
results.OLS.mixed.designmatrixPieces.constant = cat(2,zeros(size(condDesign)),constantDesign);
% close all
% figure('WindowStyle','docked'); colormap gray
% imagesc(results.OLS.mixed.designmatrix)
% figure('WindowStyle','docked'); colormap gray
% imagesc(results.OLS.mixed.designmatrixPieces.cond)
% figure('WindowStyle','docked'); colormap gray
% imagesc(results.OLS.mixed.designmatrixPieces.motion)
% figure('WindowStyle','docked'); colormap gray
% imagesc(results.OLS.mixed.designmatrixPieces.constant,[-1 1])

%% Estimate parameters with OLS
display('computing OLS')
results.OLS.mixed.parameters = ...
    mtimescell(olsmatrix2(results.OLS.mixed.designmatrix), ...
    cellfun(@(x) squish(x,dimdata)',data,'UniformOutput',0));  % regressors x voxels

designMatrixReduced = results.OLS.mixed.designmatrix;
designMatrixReduced(:,1:numruns*2) = zeros(size(designMatrixReduced(:,1:numruns*2)));

% Compute predictions
% fullModel = zeros(numvoxels,size(designMatrixFull,1));
reducedModel = zeros(numvoxels,size(designMatrixReduced,1));
for vox = 1:numvoxels
    tmpModel = bsxfun(@times,results.OLS.mixed.parameters(:,vox)',designMatrixReduced);
    reducedModel(vox,:) = sum(tmpModel,2);
end

reducedModel = reshape(reducedModel,  [xyzsize size(reducedModel,2)]);
results.OLS.mixed.parameters = reshape(results.OLS.mixed.parameters',  [xyzsize size(results.OLS.mixed.parameters',2)]);




%Subtract trend and movement from data
for runInd = 1:numruns
    tmp{runInd} = reducedModel(:,:,:,numtime*(runInd-1)+1:numtime*(runInd));
end
reducedModel = tmp;

for runInd = 1:numruns
    detrendData{runInd} = data{runInd} - reducedModel{runInd};
end



%% Convert to percent BOLD
results.OLS.mixed.constant.brain = results.OLS.mixed.parameters(:,:,:,any(results.OLS.mixed.designmatrixPieces.constant,1)).*(vectorlength(ones(numtime,1))/numtime); %vectorlength(ones(numtime,1))/numtime is the value in the design matrix
results.OLS.mixed.runInd = reshape(repmat(1:size(results.OLS.mixed.constant.brain,4),[numruns*2/size(results.OLS.mixed.constant.brain,4) 1]),[1 numruns*2]);
for run = 1:size(results.OLS.mixed.constant.brain,4)
    con = 1./abs(results.OLS.mixed.constant.brain(:,:,:,run)) * 100;
    curInd = find(results.OLS.mixed.runInd==run);
    results.OLS.mixed.parameters(:,:,:,curInd) = bsxfun(@times,results.OLS.mixed.parameters(:,:,:,curInd),con);
    detrendData{runInd} = bsxfun(@times,detrendData{runInd},con);
end

%% Compute amp and delay
[results.OLS.mixed.delay, results.OLS.mixed.amp] = cart2pol(results.OLS.mixed.parameters(:,:,:,2:2:numruns*2),results.OLS.mixed.parameters(:,:,:,1:2:numruns*2));






