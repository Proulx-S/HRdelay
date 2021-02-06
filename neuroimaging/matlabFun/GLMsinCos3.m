function [results,ResidReduced] = GLMsinCos3(design,data,stimdur,tr,hrfmodel,hrfknobs,opt,splitIn)

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





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MIXED(random)-EFFECT FIT
tic
fprintf('*** MIXED(random)-EFFECT FIT ***\n');
%% Construct design matrix
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
sessLabel = [opt.sessionLabel{:}]';

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
%deal with higher order poly
ind = find(any(results.OLS.mixed.designmatrixPieces.constant,1));
for i = 1:unique(opt.maxpolydeg)
    results.OLS.mixed.designmatrixPieces.(['p' num2str(i)]) = zeros(size(results.OLS.mixed.designmatrixPieces.constant));
    results.OLS.mixed.designmatrixPieces.(['p' num2str(i)])(:,ind+i) = results.OLS.mixed.designmatrix(:,ind+i);
end

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

%% Compute F
display('computing F stats')
% Set-up full and reduced design matrices
designMatrixFull = results.OLS.mixed.designmatrix;
designMatrixReduced = designMatrixFull;
% designMatrixReduced(:,numruns*2/3*2+1:numruns*2) = zeros(size(designMatrixReduced(:,numruns*2/3*2+1:numruns*2)));
designMatrixReduced(:,1:numruns*2) = zeros(size(designMatrixReduced(:,1:numruns*2)));
results.OLS.mixed.designMatrixReduced = designMatrixReduced;

% Compute predictions
fullModel = zeros(numvoxels,size(designMatrixFull,1));
reducedModel = zeros(numvoxels,size(designMatrixReduced,1));
for vox = 1:numvoxels
    tmpModel = bsxfun(@times,results.OLS.mixed.parameters(:,vox)',designMatrixFull);
    fullModel(vox,:) = sum(tmpModel,2);
    reducedModel(vox,:) = sum(tmpModel(:,any(designMatrixReduced,1)),2);
end
fullModel = reshape(fullModel,  [xyzsize size(fullModel,2)]);
reducedModel = reshape(reducedModel,  [xyzsize size(reducedModel,2)]);

% Compute RSSs
RSSfull = bsxfun(@minus, fullModel, catcell(4,data));
ResidFull = squeeze(mat2cell(RSSfull,xyzsize(1),xyzsize(2),xyzsize(3),ones(numruns,1).*numtime))';
RSSfull = bsxfun(@times, RSSfull, RSSfull);
RSSfull_cond3 = sum(RSSfull(:,:,:,condLabel==3),4);
RSSfull_sess1 = sum(RSSfull(:,:,:,sessLabel==1),4);
RSSfull_sess2 = sum(RSSfull(:,:,:,sessLabel==2),4);
RSSfull_cond3sess1 = sum(RSSfull(:,:,:,condLabel==3&sessLabel==1),4);
RSSfull_cond3sess2 = sum(RSSfull(:,:,:,condLabel==3&sessLabel==2),4);
RSSfull = sum(RSSfull,4);

RSSreduced = bsxfun(@minus, reducedModel, catcell(4,data));
ResidReduced = squeeze(mat2cell(RSSreduced,xyzsize(1),xyzsize(2),xyzsize(3),ones(numruns,1).*numtime))';
RSSreduced = bsxfun(@times, RSSreduced, RSSreduced);
RSSreduced_cond3 = sum(RSSreduced(:,:,:,condLabel==3),4);
RSSreduced_sess1 = sum(RSSreduced(:,:,:,sessLabel==1),4);
RSSreduced_sess2 = sum(RSSreduced(:,:,:,sessLabel==2),4);
RSSreduced_cond3sess1 = sum(RSSreduced(:,:,:,condLabel==3&sessLabel==1),4);
RSSreduced_cond3sess2 = sum(RSSreduced(:,:,:,condLabel==3&sessLabel==2),4);
RSSreduced = sum(RSSreduced,4);


% Compute and output F (considering all data)
%Define degrees of freedom
pFull = length(find(any(designMatrixFull,1)));
pReduced = pFull-2*(numruns/3);
n = size(designMatrixFull,1);
%Compute F
results.OLS.mixed.F.val.F = ((RSSreduced-RSSfull)./(pFull-pReduced)) ./ (RSSfull./(n-pFull));
results.OLS.mixed.F.val.F(results.OLS.mixed.F.val.F<=0) = eps(class(results.OLS.mixed.F.val.F)); % replace with smallest possible number if impossible value
results.OLS.mixed.F.df.pFull = pFull;
results.OLS.mixed.F.df.pReduced = pReduced;
results.OLS.mixed.F.df.n = n;
results.OLS.mixed.F.info = 'considering all data';


% Compute and output F (considering sess1 only)
%Define degrees of freedom
pFull = length(find(any(designMatrixFull(sessLabel==1,:),1)));
pReduced = pFull-2*(numruns/3/2);
n = size(designMatrixFull(sessLabel==1,:),1);
%Compute F
results.OLS.mixed.Fsess1.val.F = ((RSSreduced_sess1-RSSfull_sess1)./(pFull-pReduced)) ./ (RSSfull_sess1./(n-pFull));
results.OLS.mixed.Fsess1.val.F(results.OLS.mixed.Fsess1.val.F<=0) = eps(class(results.OLS.mixed.Fsess1.val.F)); % replace with smallest possible number if impossible value
results.OLS.mixed.Fsess1.df.pFull = pFull;
results.OLS.mixed.Fsess1.df.pReduced = pReduced;
results.OLS.mixed.Fsess1.df.n = n;
results.OLS.mixed.Fsess1.info = 'considering sess1 only';
results.OLS.mixed.Fsess1_DMind = sessLabel==1;


% Compute and output F (considering sess2 only)
%Define degrees of freedom
pFull = length(find(any(designMatrixFull(sessLabel==2,:),1)));
pReduced = pFull-2*(numruns/3/2);
n = size(designMatrixFull(sessLabel==2,:),1);
%Compute F
results.OLS.mixed.Fsess2.val.F = ((RSSreduced_sess2-RSSfull_sess2)./(pFull-pReduced)) ./ (RSSfull_sess2./(n-pFull));
results.OLS.mixed.Fsess2.val.F(results.OLS.mixed.Fsess2.val.F<=0) = eps(class(results.OLS.mixed.Fsess2.val.F)); % replace with smallest possible number if impossible value
results.OLS.mixed.Fsess2.df.pFull = pFull;
results.OLS.mixed.Fsess2.df.pReduced = pReduced;
results.OLS.mixed.Fsess2.df.n = n;
results.OLS.mixed.Fsess2.info = 'considering sess2 only';
results.OLS.mixed.Fsess2_DMind = sessLabel==2;


% Compute and output F (considering only cond3)
%Define degrees of freedom
pFull = length(find(any(designMatrixFull(condLabel==3,:),1)));
pReduced = pFull-2*(numruns/3);
n = size(designMatrixFull(condLabel==3,:),1);
%Compute F
results.OLS.mixed.Fcond3.val.F = ((RSSreduced_cond3-RSSfull_cond3)./(pFull-pReduced)) ./ (RSSfull_cond3./(n-pFull));
results.OLS.mixed.Fcond3.val.F(results.OLS.mixed.Fcond3.val.F<=0) = eps(class(results.OLS.mixed.Fcond3.val.F)); % replace with smallest possible number if impossible value
results.OLS.mixed.Fcond3.df.pFull = pFull;
results.OLS.mixed.Fcond3.df.pReduced = pReduced;
results.OLS.mixed.Fcond3.df.n = n;
results.OLS.mixed.Fcond3_DMind = condLabel==3;


% Compute and output F (considering only cond3 sess1)
%Define degrees of freedom
pFull = length(find(any(designMatrixFull(condLabel==3&sessLabel==1,:),1)));
pReduced = pFull-2*(numruns/3/2);
n = size(designMatrixFull(condLabel==3&sessLabel==1,:),1);
%Compute F
results.OLS.mixed.Fcond3sess1.val.F = ((RSSreduced_cond3sess1-RSSfull_cond3sess1)./(pFull-pReduced)) ./ (RSSfull_cond3sess1./(n-pFull));
results.OLS.mixed.Fcond3sess1.val.F(results.OLS.mixed.Fcond3sess1.val.F<=0) = eps(class(results.OLS.mixed.Fcond3sess1.val.F)); % replace with smallest possible number if impossible value
results.OLS.mixed.Fcond3sess1.df.pFull = pFull;
results.OLS.mixed.Fcond3sess1.df.pReduced = pReduced;
results.OLS.mixed.Fcond3sess1.df.n = n;
results.OLS.mixed.Fcond3sess1_DMind = condLabel==3&sessLabel==1;


% Compute and output F (considering only cond3 sess2)
%Define degrees of freedom
pFull = length(find(any(designMatrixFull(condLabel==3&sessLabel==2,:),1)));
pReduced = pFull-2*(numruns/3/2);
n = size(designMatrixFull(condLabel==3&sessLabel==2,:),1);
%Compute F
results.OLS.mixed.Fcond3sess2.val.F = ((RSSreduced_cond3sess2-RSSfull_cond3sess2)./(pFull-pReduced)) ./ (RSSfull_cond3sess2./(n-pFull));
results.OLS.mixed.Fcond3sess2.val.F(results.OLS.mixed.Fcond3sess2.val.F<=0) = eps(class(results.OLS.mixed.Fcond3sess2.val.F)); % replace with smallest possible number if impossible value
results.OLS.mixed.Fcond3sess2.df.pFull = pFull;
results.OLS.mixed.Fcond3sess2.df.pReduced = pReduced;
results.OLS.mixed.Fcond3sess2.df.n = n;
results.OLS.mixed.Fcond3sess2_DMind = condLabel==3&sessLabel==2;


clear temp temp2 sumsq good X fullModel reducedModel RSSfull RSSreduced

%% Reformat outputs
results.OLS.mixed.parameters =   reshape(results.OLS.mixed.parameters',  [xyzsize size(results.OLS.mixed.parameters,1)]);
results.OLS.mixed.poly.p0 = results.OLS.mixed.parameters(:,:,:,any(results.OLS.mixed.designmatrixPieces.constant,1));
for i = 1:unique(opt.maxpolydeg)
    results.OLS.mixed.poly.(['p' num2str(i)]) = results.OLS.mixed.parameters(:,:,:,any(results.OLS.mixed.designmatrixPieces.(['p' num2str(i)]),1));
end
results.OLS.mixed.runInd = reshape(repmat(1:size(results.OLS.mixed.poly.p0,4),[numruns*2/size(results.OLS.mixed.poly.p0,4) 1]),[1 numruns*2]);
% for run = 1:size(results.OLS.mixed.poly.p0,4)
%     con = 1./abs(results.OLS.mixed.poly.p0(:,:,:,run)) * 100;
%     curInd = find(results.OLS.mixed.runInd==run);
%     results.OLS.mixed.parameters(:,:,:,curInd) = bsxfun(@times,results.OLS.mixed.parameters(:,:,:,curInd),con);
% end

%% Compute amp and delay
[results.OLS.mixed.delay, results.OLS.mixed.amp] = cart2pol(results.OLS.mixed.parameters(:,:,:,1:2:numruns*2),results.OLS.mixed.parameters(:,:,:,2:2:numruns*2));
toc

tic
% %% Directly convert the data to %BOLD
% % Now that we have mean BOLD properly estimated, doing this conversion
% % directly on the data will avoid some problems in the fixed effect model,
% % where conversion to %BOLD on the parameter estimate is suboptimal since
% % each parameter is estimated on data with heterogeneous mean BOLD.
% fprintf('*** CONVERTING BOLD TIME SERIES TO PERCENT BOLD ***\n');
% for run = 1:size(results.OLS.mixed.poly.p0,4)
%     oneOver_curBrain = 1./results.OLS.mixed.poly.p0(:,:,:,run);
%     
%     curData = data{results.OLS.mixed.runInd(1:2:end)==run};
%     curData = bsxfun(@times,curData,oneOver_curBrain);
%     data{results.OLS.mixed.runInd(1:2:end)==run} = curData;
% end

fprintf('*** COMPUTE VEIN (SD of residual / signal baseline) ***\n');

results.OLS.mixed.veinFull = nan([xyzsize numruns]);
results.OLS.mixed.veinFull_info = 'std of full model (including regressors of interest) residuals  /  constant term';
results.OLS.mixed.veinReduced = nan([xyzsize numruns]);
results.OLS.mixed.veinReduced_info = 'std of reduced model (excluding regressors of interest) residuals  /  constant term';
for run = 1:size(results.OLS.mixed.poly.p0,4)
    base = results.OLS.mixed.poly.p0(:,:,:,run);
    
    sd = std(ResidFull{run},[],4);
    results.OLS.mixed.veinFull(:,:,:,run) = sd./base;
    
    sd = std(ResidReduced{run},[],4);
    results.OLS.mixed.veinReduced(:,:,:,run) = sd./base;
end

toc



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




