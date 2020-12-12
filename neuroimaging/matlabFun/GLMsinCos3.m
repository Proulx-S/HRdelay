function [results,sessModel] = GLMsinCos3(design,data,stimdur,tr,hrfmodel,hrfknobs,opt,splitIn)

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





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MIXED-EFFECT FIT
tic
fprintf('*** MIXED-EFFECT FIT ***\n');
%% Construct design matrix
poly = cell(1,numruns);
polyNmotion = cell(1,numruns);
constant = cell(1,numruns);

numtime = size(data{p},dimtime)*splitIn;
for p=1:numruns/splitIn
    tmp = constructpolynomialmatrix(numtime,0:opt.maxpolydeg(p));
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

%% Compute F
display('computing F stats')
% Set-up full and reduced design matrices
designMatrixFull = results.OLS.mixed.designmatrix;
designMatrixReduced = designMatrixFull;
designMatrixReduced(:,numruns*2/3*2+1:numruns*2) = zeros(size(designMatrixReduced(:,numruns*2/3*2+1:numruns*2)));
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

% Find indices of data to include for different F tests
ind_cond3 = size(fullModel,4)*2/3+1:size(fullModel,4);
ind_sess1 = [size(fullModel,4)*0/6+1:size(fullModel,4)*1/6 size(fullModel,4)*2/6+1:size(fullModel,4)*3/6 size(fullModel,4)*4/6+1:size(fullModel,4)*5/6];
ind_sess2 = [size(fullModel,4)*1/6+1:size(fullModel,4)*2/6 size(fullModel,4)*3/6+1:size(fullModel,4)*4/6 size(fullModel,4)*5/6+1:size(fullModel,4)*6/6];
ind_cond3sess1 = size(fullModel,4)*4/6+1:size(fullModel,4)*5/6;
ind_cond3sess2 = size(fullModel,4)*5/6+1:size(fullModel,4)*6/6;
% figure('WindowStyle','docked');
% tmp = nan(size(designMatrixFull,1),1);
% tmp(ind_sess1) = 1;
% imagesc(cat(2,tmp,designMatrixFull))

% Compute RSSs
RSSfull = bsxfun(@minus, fullModel, catcell(4,data));
ResidFull = squeeze(mat2cell(RSSfull,xyzsize(1),xyzsize(2),xyzsize(3),ones(numruns,1).*numtime))';
RSSfull = bsxfun(@times, RSSfull, RSSfull);
RSSfull_cond3 = sum(RSSfull(:,:,:,ind_cond3),4);
RSSfull_sess1 = sum(RSSfull(:,:,:,ind_sess1),4);
RSSfull_sess2 = sum(RSSfull(:,:,:,ind_sess2),4);
RSSfull_cond3sess1 = sum(RSSfull(:,:,:,ind_cond3sess1),4);
RSSfull_cond3sess2 = sum(RSSfull(:,:,:,ind_cond3sess2),4);
RSSfull = sum(RSSfull,4);

RSSreduced = bsxfun(@minus, reducedModel, catcell(4,data));
% ResidReduced = squeeze(mat2cell(RSSreduced,xyzsize(1),xyzsize(2),xyzsize(3),ones(numruns,1).*numtime))';
RSSreduced = bsxfun(@times, RSSreduced, RSSreduced);
RSSreduced_cond3 = sum(RSSreduced(:,:,:,ind_cond3),4);
RSSreduced_sess1 = sum(RSSreduced(:,:,:,ind_sess1),4);
RSSreduced_sess2 = sum(RSSreduced(:,:,:,ind_sess2),4);
RSSreduced_cond3sess1 = sum(RSSreduced(:,:,:,ind_cond3sess1),4);
RSSreduced_cond3sess2 = sum(RSSreduced(:,:,:,ind_cond3sess2),4);
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
pFull = length(find(any(designMatrixFull(ind_sess1,:),1)));
pReduced = pFull-2*(numruns/3/2);
n = size(designMatrixFull(ind_sess1,:),1);
%Compute F
results.OLS.mixed.Fsess1.val.F = ((RSSreduced_sess1-RSSfull_sess1)./(pFull-pReduced)) ./ (RSSfull_sess1./(n-pFull));
results.OLS.mixed.Fsess1.val.F(results.OLS.mixed.Fsess1.val.F<=0) = eps(class(results.OLS.mixed.Fsess1.val.F)); % replace with smallest possible number if impossible value
results.OLS.mixed.Fsess1.df.pFull = pFull;
results.OLS.mixed.Fsess1.df.pReduced = pReduced;
results.OLS.mixed.Fsess1.df.n = n;
results.OLS.mixed.Fsess1.info = 'considering sess1 only';
results.OLS.mixed.Fsess1_DMind = ind_sess1;


% Compute and output F (considering sess2 only)
%Define degrees of freedom
pFull = length(find(any(designMatrixFull(ind_sess2,:),1)));
pReduced = pFull-2*(numruns/3/2);
n = size(designMatrixFull(ind_sess2,:),1);
%Compute F
results.OLS.mixed.Fsess2.val.F = ((RSSreduced_sess2-RSSfull_sess2)./(pFull-pReduced)) ./ (RSSfull_sess2./(n-pFull));
results.OLS.mixed.Fsess2.val.F(results.OLS.mixed.Fsess2.val.F<=0) = eps(class(results.OLS.mixed.Fsess2.val.F)); % replace with smallest possible number if impossible value
results.OLS.mixed.Fsess2.df.pFull = pFull;
results.OLS.mixed.Fsess2.df.pReduced = pReduced;
results.OLS.mixed.Fsess2.df.n = n;
results.OLS.mixed.Fsess2.info = 'considering sess2 only';
results.OLS.mixed.Fsess2_DMind = ind_sess2;


% Compute and output F (considering only cond3)
%Define degrees of freedom
pFull = length(find(any(designMatrixFull(ind_cond3,:),1)));
pReduced = pFull-2*(numruns/3);
n = size(designMatrixFull(ind_cond3,:),1);
%Compute F
results.OLS.mixed.Fcond3.val.F = ((RSSreduced_cond3-RSSfull_cond3)./(pFull-pReduced)) ./ (RSSfull_cond3./(n-pFull));
results.OLS.mixed.Fcond3.val.F(results.OLS.mixed.Fcond3.val.F<=0) = eps(class(results.OLS.mixed.Fcond3.val.F)); % replace with smallest possible number if impossible value
results.OLS.mixed.Fcond3.df.pFull = pFull;
results.OLS.mixed.Fcond3.df.pReduced = pReduced;
results.OLS.mixed.Fcond3.df.n = n;
results.OLS.mixed.Fcond3_DMind = ind_cond3;


% Compute and output F (considering only cond3 sess1)
%Define degrees of freedom
pFull = length(find(any(designMatrixFull(ind_cond3sess1,:),1)));
pReduced = pFull-2*(numruns/3/2);
n = size(designMatrixFull(ind_cond3sess1,:),1);
%Compute F
results.OLS.mixed.Fcond3sess1.val.F = ((RSSreduced_cond3sess1-RSSfull_cond3sess1)./(pFull-pReduced)) ./ (RSSfull_cond3sess1./(n-pFull));
results.OLS.mixed.Fcond3sess1.val.F(results.OLS.mixed.Fcond3sess1.val.F<=0) = eps(class(results.OLS.mixed.Fcond3sess1.val.F)); % replace with smallest possible number if impossible value
results.OLS.mixed.Fcond3sess1.df.pFull = pFull;
results.OLS.mixed.Fcond3sess1.df.pReduced = pReduced;
results.OLS.mixed.Fcond3sess1.df.n = n;
results.OLS.mixed.Fcond3sess1_DMind = ind_cond3sess1;


% Compute and output F (considering only cond3 sess2)
%Define degrees of freedom
pFull = length(find(any(designMatrixFull(ind_cond3sess2,:),1)));
pReduced = pFull-2*(numruns/3/2);
n = size(designMatrixFull(ind_cond3sess2,:),1);
%Compute F
results.OLS.mixed.Fcond3sess2.val.F = ((RSSreduced_cond3sess2-RSSfull_cond3sess2)./(pFull-pReduced)) ./ (RSSfull_cond3sess2./(n-pFull));
results.OLS.mixed.Fcond3sess2.val.F(results.OLS.mixed.Fcond3sess2.val.F<=0) = eps(class(results.OLS.mixed.Fcond3sess2.val.F)); % replace with smallest possible number if impossible value
results.OLS.mixed.Fcond3sess2.df.pFull = pFull;
results.OLS.mixed.Fcond3sess2.df.pReduced = pReduced;
results.OLS.mixed.Fcond3sess2.df.n = n;
results.OLS.mixed.Fcond3sess2_DMind = ind_cond3sess2;


clear temp temp2 sumsq good X fullModel reducedModel RSSfull RSSreduced

%% Reformat outputs
results.OLS.mixed.parameters =   reshape(results.OLS.mixed.parameters',  [xyzsize size(results.OLS.mixed.parameters,1)]);

%% Convert to percent BOLD
results.OLS.mixed.constant.brain = results.OLS.mixed.parameters(:,:,:,any(results.OLS.mixed.designmatrixPieces.constant,1)).*(vectorlength(ones(numtime,1))/numtime); %vectorlength(ones(numtime,1))/numtime is the value in the design matrix
results.OLS.mixed.runInd = reshape(repmat(1:size(results.OLS.mixed.constant.brain,4),[numruns*2/size(results.OLS.mixed.constant.brain,4) 1]),[1 numruns*2]);
for run = 1:size(results.OLS.mixed.constant.brain,4)
    con = 1./abs(results.OLS.mixed.constant.brain(:,:,:,run)) * 100;
    curInd = find(results.OLS.mixed.runInd==run);
    results.OLS.mixed.parameters(:,:,:,curInd) = bsxfun(@times,results.OLS.mixed.parameters(:,:,:,curInd),con);
end

%% Compute amp and delay
[results.OLS.mixed.delay, results.OLS.mixed.amp] = cart2pol(results.OLS.mixed.parameters(:,:,:,1:2:numruns*2),results.OLS.mixed.parameters(:,:,:,2:2:numruns*2));
toc

tic
%% Directly convert the data to %BOLD
% Now that we have mean BOLD properly estimated, doing this conversion
% directly on the data will avoid some problems in the fixed effect model,
% where conversion to %BOLD on the parameter estimate is suboptimal since
% each parameter is estimated on data with heterogeneous mean BOLD.
fprintf('*** CONVERTING BOLD TIME SERIES TO PERCENT BOLD ***\n');
for run = 1:size(results.OLS.mixed.constant.brain,4)
    oneOver_curBrain = 1./results.OLS.mixed.constant.brain(:,:,:,run);
    
    curData = data{results.OLS.mixed.runInd(1:2:end)==run};
    curData = bsxfun(@times,curData,oneOver_curBrain);
    data{results.OLS.mixed.runInd(1:2:end)==run} = curData;
end

fprintf('*** COMPUTE VEIN (SD of residual / signal baseline) ***\n');

results.OLS.mixed.vein = nan([xyzsize numruns]);
for run = 1:size(results.OLS.mixed.constant.brain,4)
    sd = std(ResidFull{run},[],4);
    base = results.OLS.mixed.constant.brain(:,:,:,run);
    results.OLS.mixed.vein(:,:,:,run) = sd./base;
end

toc


if opt.splitedIn==1 || ~isfield(opt,'sessModel')
%     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIXED-EFFECT FIT
%     tic
%     fprintf('*** FIXED-EFFECT FIT ***\n');
%     %% Construct design matrix
%     % Polynomials and motion
%     poly = cell(1,numruns);
%     polyNmotion = cell(1,numruns);
%     constant = cell(1,numruns);
%     
%     numtime = size(data{p},dimtime)*splitIn;
%     for p=1:numruns/splitIn
%         tmp = constructpolynomialmatrix(numtime,0:opt.maxpolydeg(p));
%         poly((p-1)*splitIn+1:p*splitIn) = mat2cell(tmp,repmat(size(data{p},dimtime),1,splitIn),2)';
%     end
%     
%     numtime = size(data{p},dimtime)*splitIn;
%     for p=1:numruns
%         polyNmotion{p} = cat(2,poly{p}, ...
%             opt.extraregressors{p}, ...
%             []);
%         tmpConstant = poly{p}; tmpConstant(:,2:end) = zeros(size(tmpConstant(:,2:end)));
%         constant{p} = cat(2,tmpConstant, ...
%             zeros(size(opt.extraregressors{p})), ...
%             []);
%     end
%     %Reformat to remerge motion regressors that were split before
%     polyNmotion_b = reshape(polyNmotion,[opt.splitedIn length(polyNmotion)/opt.splitedIn])';
%     constant_b = reshape(constant,[opt.splitedIn length(constant)/opt.splitedIn])';
%     polyNmotion = cell(1,size(polyNmotion_b,1));
%     constant = cell(1,size(polyNmotion_b,1));
%     for ii = 1:size(polyNmotion_b,1)
%         polyNmotion{ii} = catcell(1,polyNmotion_b(ii,:));
%         constant{ii} = catcell(1,constant_b(ii,:));
%     end
%     clear polyNmotion_b constant_b
%     
%     % Conditions
%     condLabel = cat(1,ones(numruns/3,1)*1,ones(numruns/3,1)*2,ones(numruns/3,1)*3);
%     
%     condDesign = cell(1,3);
%     polyNmotionDesign = cell(1,3);
%     constantDesign = cell(1,3);
%     for cond = 1:3
%         condInd = condLabel==cond;
%         condInd_polyNmotion = condInd(1:opt.splitedIn:end);
%         condDesign{cond} = catcell(1,cache.rawdesign(condInd));
%         
%         polyNmotionDesign{cond} = blkdiag(polyNmotion{condInd_polyNmotion});
%         constantDesign{cond} = blkdiag(constant{condInd_polyNmotion});
%     end
%     
%     % Merge conditions, polynomials and motion
%     condDesign = blkdiag(condDesign{:});
%     % sessDesign = catcell(1,sessDesign); sessDesign(:,1:2) = [];
%     polyNmotionDesign = blkdiag(polyNmotionDesign{:});
%     constantDesign = blkdiag(constantDesign{:});
%     
%     results.OLS.fixed.designmatrix = cat(2,condDesign,polyNmotionDesign);% time x regressors
%     % results.OLS.fixed.designmatrixCond3 = cat(2,condDesign(end*2/3+1:end,:),polyNmotionDesign(end*2/3+1:end,:));% time x regressors
%     
%     % Extract some usefull params
%     results.OLS.fixed.designmatrixPieces.cond = cat(2,condDesign,zeros(size(polyNmotionDesign)));
%     results.OLS.fixed.designmatrixPieces.motion = cat(2,zeros(size(condDesign)),polyNmotionDesign);
%     results.OLS.fixed.designmatrixPieces.constant = cat(2,zeros(size(condDesign)),constantDesign);
% %     close all
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.fixed.designmatrix)
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.fixed.designmatrixPieces.cond)
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.fixed.designmatrixPieces.motion)
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.fixed.designmatrixPieces.constant,[-1 1])
%     
%     
%     %% Estimate parameters with OLS
%     display('computing OLS')
%     results.OLS.fixed.parameters = ...
%         mtimescell(olsmatrix2(results.OLS.fixed.designmatrix), ...
%         cellfun(@(x) squish(x,dimdata)',data,'UniformOutput',0));  % regressors x voxels
%     
%     %% Compute F for cond 3
%     display('computing F stats for cond 3')
%     % Set-up full and reduced design matrices
%     designMatrixFull = results.OLS.fixed.designmatrix;
%     designMatrixReduced = designMatrixFull;
%     designMatrixReduced(:,5:6) = zeros(size(designMatrixReduced(:,5:6)));
%     % close all
%     % figure('WindowStyle','docked'); colormap gray
%     % imagesc(designMatrixFull)
%     % figure('WindowStyle','docked'); colormap gray
%     % imagesc(designMatrixReduced)
%     
%     % Compute predictions
%     fullModel = zeros(numvoxels,size(designMatrixFull,1));
%     reducedModel = zeros(numvoxels,size(designMatrixReduced,1));
%     for vox = 1:numvoxels
%         tmpModel = bsxfun(@times,results.OLS.fixed.parameters(:,vox)',designMatrixFull);
%         fullModel(vox,:) = sum(tmpModel,2);
%         reducedModel(vox,:) = sum(tmpModel(:,any(designMatrixReduced,1)),2);
%     end
%     fullModel = reshape(fullModel,  [xyzsize size(fullModel,2)]);
%     reducedModel = reshape(reducedModel,  [xyzsize size(reducedModel,2)]);
%     
%     
%     % Find indices of data to include for different F tests
%     ind_cond3 = size(fullModel,4)*2/3+1:size(fullModel,4);
%     ind_sess1 = [size(fullModel,4)*0/6+1:size(fullModel,4)*1/6 size(fullModel,4)*2/6+1:size(fullModel,4)*3/6 size(fullModel,4)*4/6+1:size(fullModel,4)*5/6];
%     ind_sess2 = [size(fullModel,4)*1/6+1:size(fullModel,4)*2/6 size(fullModel,4)*3/6+1:size(fullModel,4)*4/6 size(fullModel,4)*5/6+1:size(fullModel,4)*6/6];
%     ind_cond3sess1 = size(fullModel,4)*4/6+1:size(fullModel,4)*5/6;
%     ind_cond3sess2 = size(fullModel,4)*5/6+1:size(fullModel,4)*6/6;
%     
% %     figure('WindowStyle','docked'); colormap gray
% %     tmp = zeros(size(fullModel,4),1);
% %     tmp(ind_cond3sess2) = 1;
% %     imagesc(cat(2,tmp,results.OLS.fixed.designmatrix)); xlim([0.5 50])
%     
%     % Compute RSSs
%     RSSfull = bsxfun(@minus, fullModel, catcell(4,data));
%     RSSfull = bsxfun(@times, RSSfull, RSSfull);
%     RSSfull_cond3 = sum(RSSfull(:,:,:,ind_cond3),4);
%     RSSfull_sess1 = sum(RSSfull(:,:,:,ind_sess1),4);
%     RSSfull_sess2 = sum(RSSfull(:,:,:,ind_sess2),4);
%     RSSfull_cond3sess1 = sum(RSSfull(:,:,:,ind_cond3sess1),4);
%     RSSfull_cond3sess2 = sum(RSSfull(:,:,:,ind_cond3sess2),4);
%     RSSfull = sum(RSSfull,4);
%     
%     RSSreduced = bsxfun(@minus, reducedModel, catcell(4,data));
%     RSSreduced = bsxfun(@times, RSSreduced, RSSreduced);
%     RSSreduced_cond3 = sum(RSSreduced(:,:,:,ind_cond3),4);
%     RSSreduced_sess1 = sum(RSSreduced(:,:,:,ind_sess1),4);
%     RSSreduced_sess2 = sum(RSSreduced(:,:,:,ind_sess2),4);
%     RSSreduced_cond3sess1 = sum(RSSreduced(:,:,:,ind_cond3sess1),4);
%     RSSreduced_cond3sess2 = sum(RSSreduced(:,:,:,ind_cond3sess2),4);
%     RSSreduced = sum(RSSreduced,4);
%     
%     % Compute and output F (considering all data)
%     %Define degrees of freedom
%     pFull = length(find(any(designMatrixFull,1)));
%     pReduced = pFull-2;
%     n = size(designMatrixFull,1);
%     %Compute F
%     results.OLS.fixed.F.val.F = ((RSSreduced-RSSfull)./(pFull-pReduced)) ./ (RSSfull./(n-pFull));
%     results.OLS.fixed.F.val.F(results.OLS.fixed.F.val.F<=0) = eps(class(results.OLS.fixed.F.val.F)); % replace with smallest possible number if impossible value
%     results.OLS.fixed.F.df.pFull = pFull;
%     results.OLS.fixed.F.df.pReduced = pReduced;
%     results.OLS.fixed.F.df.n = n;
%     
%     % Compute and output F (considering sess1 only)
%     %Define degrees of freedom
%     pFull = length(find(any(designMatrixFull(ind_sess1,:),1)));
%     pReduced = pFull-2;
%     n = size(designMatrixFull(ind_sess1,:),1);
%     %Compute F
%     results.OLS.fixed.Fsess1.val.F = ((RSSreduced_sess1-RSSfull_sess1)./(pFull-pReduced)) ./ (RSSfull_sess1./(n-pFull));
%     results.OLS.fixed.Fsess1.val.F(results.OLS.fixed.Fsess1.val.F<=0) = eps(class(results.OLS.fixed.Fsess1.val.F)); % replace with smallest possible number if impossible value
%     results.OLS.fixed.Fsess1.df.pFull = pFull;
%     results.OLS.fixed.Fsess1.df.pReduced = pReduced;
%     results.OLS.fixed.Fsess1.df.n = n;
%     
%     % Compute and output F (considering sess2 only)
%     %Define degrees of freedom
%     pFull = length(find(any(designMatrixFull(ind_sess2,:),1)));
%     pReduced = pFull-2;
%     n = size(designMatrixFull(ind_sess2,:),1);
%     %Compute F
%     results.OLS.fixed.Fsess2.val.F = ((RSSreduced_sess2-RSSfull_sess2)./(pFull-pReduced)) ./ (RSSfull_sess2./(n-pFull));
%     results.OLS.fixed.Fsess2.val.F(results.OLS.fixed.Fsess2.val.F<=0) = eps(class(results.OLS.fixed.Fsess2.val.F)); % replace with smallest possible number if impossible value
%     results.OLS.fixed.Fsess2.df.pFull = pFull;
%     results.OLS.fixed.Fsess2.df.pReduced = pReduced;
%     results.OLS.fixed.Fsess2.df.n = n;
%     
%     % Compute and output F (considering only cond3)
%     %Define degrees of freedom
%     pFull = length(find(any(designMatrixFull(ind_cond3,:),1)));
%     pReduced = pFull-2;
%     n = size(designMatrixFull(ind_cond3,:),1);
%     %Compute F
%     results.OLS.fixed.Fcond3.val.F = ((RSSreduced_cond3-RSSfull_cond3)./(pFull-pReduced)) ./ (RSSfull_cond3./(n-pFull));
%     results.OLS.fixed.Fcond3.val.F(results.OLS.fixed.Fcond3.val.F<=0) = eps(class(results.OLS.fixed.Fcond3.val.F)); % replace with smallest possible number if impossible value
%     results.OLS.fixed.Fcond3.df.pFull = pFull;
%     results.OLS.fixed.Fcond3.df.pReduced = pReduced;
%     results.OLS.fixed.Fcond3.df.n = n;
%     
%     % Compute and output F (considering only cond3 sess1)
%     %Define degrees of freedom
%     pFull = length(find(any(designMatrixFull(ind_cond3sess1,:),1)));
%     pReduced = pFull-2;
%     n = size(designMatrixFull(ind_cond3sess1,:),1);
%     %Compute F
%     results.OLS.fixed.Fcond3sess1.val.F = ((RSSreduced_cond3sess1-RSSfull_cond3sess1)./(pFull-pReduced)) ./ (RSSfull_cond3sess1./(n-pFull));
%     results.OLS.fixed.Fcond3sess1.val.F(results.OLS.fixed.Fcond3sess1.val.F<=0) = eps(class(results.OLS.fixed.Fcond3sess1.val.F)); % replace with smallest possible number if impossible value
%     results.OLS.fixed.Fcond3sess1.df.pFull = pFull;
%     results.OLS.fixed.Fcond3sess1.df.pReduced = pReduced;
%     results.OLS.fixed.Fcond3sess1.df.n = n;
%     
%     % Compute and output F (considering only cond3 sess2)
%     %Define degrees of freedom
%     pFull = length(find(any(designMatrixFull(ind_cond3sess2,:),1)));
%     pReduced = pFull-2;
%     n = size(designMatrixFull(ind_cond3sess2,:),1);
%     %Compute F
%     results.OLS.fixed.Fcond3sess2.val.F = ((RSSreduced_cond3sess2-RSSfull_cond3sess2)./(pFull-pReduced)) ./ (RSSfull_cond3sess2./(n-pFull));
%     results.OLS.fixed.Fcond3sess2.val.F(results.OLS.fixed.Fcond3sess2.val.F<=0) = eps(class(results.OLS.fixed.Fcond3sess2.val.F)); % replace with smallest possible number if impossible value
%     results.OLS.fixed.Fcond3sess2.df.pFull = pFull;
%     results.OLS.fixed.Fcond3sess2.df.pReduced = pReduced;
%     results.OLS.fixed.Fcond3sess2.df.n = n;
%     
%     clear temp temp2 sumsq good X fullModel reducedModel RSSfull RSSreduced
%     
%     %% Reformat outputs
%     results.OLS.fixed.parameters =   reshape(results.OLS.fixed.parameters',  [xyzsize size(results.OLS.fixed.parameters,1)]);
%     
% %     %% Convert to percent BOLD
% %     results.OLS.fixed.constant.brain = mean(results.OLS.fixed.parameters(:,:,:,any(results.OLS.fixed.designmatrixPieces.constant,1)).*(vectorlength(ones(numtime,1))/numtime),4); %vectorlength(ones(numtime,1))/numtime is the value in the design matrix
% %     con = 1./abs(results.OLS.fixed.constant.brain) * 100;
% %     results.OLS.fixed.parameters = bsxfun(@times,results.OLS.fixed.parameters,con);
%     
%     %% Compute amp and delay
%     [results.OLS.fixed.delay, results.OLS.fixed.amp] = cart2pol(results.OLS.fixed.parameters(:,:,:,1:2:6),results.OLS.fixed.parameters(:,:,:,2:2:6));
%     toc
    
    
    
    
    
    
    
    
    
%     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIXED-EFFECT FIT WITH SESSION REGRESSOR
%     tic
%     fprintf('*** FIXED-EFFECT FIT WITH SESSION REGRESSOR ***\n');
%     %% Construct design matrix
%     % Polynomials and motion
%     poly = cell(1,numruns);
%     polyNmotion = cell(1,numruns);
%     constant = cell(1,numruns);
%     
%     numtime = size(data{p},dimtime)*splitIn;
%     for p=1:numruns/splitIn
%         tmp = constructpolynomialmatrix(numtime,0:opt.maxpolydeg(p));
%         poly((p-1)*splitIn+1:p*splitIn) = mat2cell(tmp,repmat(size(data{p},dimtime),1,splitIn),2)';
%     end
%     
%     numtime = size(data{p},dimtime)*splitIn;
%     for p=1:numruns
%         polyNmotion{p} = cat(2,poly{p}, ...
%             opt.extraregressors{p}, ...
%             []);
%         tmpConstant = poly{p}; tmpConstant(:,2:end) = zeros(size(tmpConstant(:,2:end)));
%         constant{p} = cat(2,tmpConstant, ...
%             zeros(size(opt.extraregressors{p})), ...
%             []);
%     end
%     %Reformat to remerge motion regressors that were split before
%     polyNmotion_b = reshape(polyNmotion,[opt.splitedIn length(polyNmotion)/opt.splitedIn])';
%     constant_b = reshape(constant,[opt.splitedIn length(constant)/opt.splitedIn])';
%     polyNmotion = cell(1,size(polyNmotion_b,1));
%     constant = cell(1,size(polyNmotion_b,1));
%     for ii = 1:size(polyNmotion_b,1)
%         polyNmotion{ii} = catcell(1,polyNmotion_b(ii,:));
%         constant{ii} = catcell(1,constant_b(ii,:));
%     end
%     clear polyNmotion_b constant_b
%     
%     % Conditions and Sessions
%     sessionLabel = catcell(1,opt.sessionLabel);
%     condLabel = cat(1,ones(numruns/3,1)*1,ones(numruns/3,1)*2,ones(numruns/3,1)*3);
%     
%     condDesign = cell(1,3);
%     sessDesign = cell(1,3);
%     polyNmotionDesign = cell(1,3);
%     constantDesign = cell(1,3);
%     for cond = 1:3
%         condInd = condLabel==cond;
%         condInd_polyNmotion = condInd(1:opt.splitedIn:end);
%         condDesign{cond} = catcell(1,cache.rawdesign(condInd));
%         
%         curDesign = cell(1,length(unique(sessionLabel)));
%         for sess = 1:length(unique(sessionLabel))
%             sessInd = sessionLabel==sess;
%             if sess==1
%                 curDesign{sess} = zeros(size(catcell(1,cache.rawdesign(condInd&sessInd))));
%             else
%                 curDesign{sess} = catcell(1,cache.rawdesign(condInd&sessInd));
%             end
%         end
%         sessDesign{cond} = blkdiag(curDesign{:});
%         
%         polyNmotionDesign{cond} = blkdiag(polyNmotion{condInd_polyNmotion});
%         constantDesign{cond} = blkdiag(constant{condInd_polyNmotion});
%     end
%     
%     % Merge conditions, sessions, polynomials and motion
%     condDesign = blkdiag(condDesign{:});
%     sessDesign = catcell(1,sessDesign); sessDesign(:,1:2) = [];
%     polyNmotionDesign = blkdiag(polyNmotionDesign{:});
%     constantDesign = blkdiag(constantDesign{:});
%     
%     results.OLS.fixed_sessReg.designmatrix = cat(2,condDesign,sessDesign,polyNmotionDesign);% time x regressors
%     
%     % Extract some usefull params
%     results.OLS.fixed_sessReg.designmatrixPieces.cond = cat(2,condDesign,zeros(size(sessDesign)),zeros(size(polyNmotionDesign)));
%     results.OLS.fixed_sessReg.designmatrixPieces.sess = cat(2,zeros(size(condDesign)),sessDesign,zeros(size(polyNmotionDesign)));
%     results.OLS.fixed_sessReg.designmatrixPieces.motion = cat(2,zeros(size(condDesign)),zeros(size(sessDesign)),polyNmotionDesign);
%     results.OLS.fixed_sessReg.designmatrixPieces.constant = cat(2,zeros(size(condDesign)),zeros(size(sessDesign)),constantDesign);
% %     close all
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.fixed_sessReg.designmatrix)
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.fixed_sessReg.designmatrixPieces.cond)
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.fixed_sessReg.designmatrixPieces.sess)
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.fixed_sessReg.designmatrixPieces.motion)
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.fixed_sessReg.designmatrixPieces.constant,[-1 1])
%     
%     %% Estimate parameters with OLS
%     display('computing OLS')
%     results.OLS.fixed_sessReg.parameters = ...
%         mtimescell(olsmatrix2(results.OLS.fixed_sessReg.designmatrix), ...
%         cellfun(@(x) squish(x,dimdata)',data,'UniformOutput',0));  % regressors x voxels
%     
%     %% Compute F for cond 3
%     display('computing F stats for cond 3')
%     
%     % Set-up full and reduced design matrices
%     designMatrixFull = results.OLS.fixed_sessReg.designmatrix;
%     designMatrixReduced = designMatrixFull;
%     designMatrixReduced(:,5:6) = zeros(size(designMatrixReduced(:,5:6)));
%     % close all
%     % figure('WindowStyle','docked'); colormap gray
%     % imagesc(designMatrixFull)
%     % figure('WindowStyle','docked'); colormap gray
%     % imagesc(designMatrixReduced)
%     
%     % Compute predictions
%     fullModel = zeros(numvoxels,size(designMatrixFull,1));
%     reducedModel = zeros(numvoxels,size(designMatrixFull,1));
%     sessModel = zeros(numvoxels,size(designMatrixFull,1));
%     for vox = 1:numvoxels
%         tmpModel = bsxfun(@times,results.OLS.fixed_sessReg.parameters(:,vox)',designMatrixFull);
%         fullModel(vox,:) = sum(tmpModel,2);
%         reducedModel(vox,:) = sum(tmpModel(:,any(designMatrixReduced,1)),2);
%         sessModel(vox,:) = sum(tmpModel(:,any(results.OLS.fixed_sessReg.designmatrixPieces.sess,1)),2);
%     end
%     fullModel = reshape(fullModel,  [xyzsize size(fullModel,2)]);
%     reducedModel = reshape(reducedModel,  [xyzsize size(reducedModel,2)]);
%     sessModel = reshape(sessModel,  [xyzsize size(sessModel,2)]);
%     sessModel = squeeze(mat2cell(sessModel,size(sessModel,1),size(sessModel,2),size(sessModel,3),repmat(size(sessModel,4)/numruns,1,numruns)))'; % for this to match format of data, since we will subtract sessModel from data
%     
%     % Find indices of data to include for different F tests
%     ind_cond3 = size(fullModel,4)*2/3+1:size(fullModel,4);
%     ind_sess1 = [size(fullModel,4)*0/6+1:size(fullModel,4)*1/6 size(fullModel,4)*2/6+1:size(fullModel,4)*3/6 size(fullModel,4)*4/6+1:size(fullModel,4)*5/6];
%     ind_sess2 = [size(fullModel,4)*1/6+1:size(fullModel,4)*2/6 size(fullModel,4)*3/6+1:size(fullModel,4)*4/6 size(fullModel,4)*5/6+1:size(fullModel,4)*6/6];
%     ind_cond3sess1 = size(fullModel,4)*4/6+1:size(fullModel,4)*5/6;
%     ind_cond3sess2 = size(fullModel,4)*5/6+1:size(fullModel,4)*6/6;
%     
%     % Compute RSSs
%     RSSfull = bsxfun(@minus, fullModel, catcell(4,data));
%     RSSfull = bsxfun(@times, RSSfull, RSSfull);
%     RSSfull_cond3 = sum(RSSfull(:,:,:,ind_cond3),4);
%     RSSfull_sess1 = sum(RSSfull(:,:,:,ind_sess1),4);
%     RSSfull_sess2 = sum(RSSfull(:,:,:,ind_sess2),4);
%     RSSfull_cond3sess1 = sum(RSSfull(:,:,:,ind_cond3sess1),4);
%     RSSfull_cond3sess2 = sum(RSSfull(:,:,:,ind_cond3sess2),4);
%     RSSfull = sum(RSSfull,4);
%     
%     RSSreduced = bsxfun(@minus, reducedModel, catcell(4,data));
%     RSSreduced = bsxfun(@times, RSSreduced, RSSreduced);
%     RSSreduced_cond3 = sum(RSSreduced(:,:,:,ind_cond3),4);
%     RSSreduced_sess1 = sum(RSSreduced(:,:,:,ind_sess1),4);
%     RSSreduced_sess2 = sum(RSSreduced(:,:,:,ind_sess2),4);
%     RSSreduced_cond3sess1 = sum(RSSreduced(:,:,:,ind_cond3sess1),4);
%     RSSreduced_cond3sess2 = sum(RSSreduced(:,:,:,ind_cond3sess2),4);
%     RSSreduced = sum(RSSreduced,4);
%     
%     % Compute and output F (considering all data)
%     %Define degrees of freedom
%     pFull = length(find(any(designMatrixFull,1)));
%     pReduced = pFull-2;
%     n = size(designMatrixFull,1);
%     %Compute F
%     results.OLS.fixed_sessReg.F.val.F = ((RSSreduced-RSSfull)./(pFull-pReduced)) ./ (RSSfull./(n-pFull));
%     results.OLS.fixed_sessReg.F.val.F(results.OLS.fixed_sessReg.F.val.F<=0) = eps(class(results.OLS.fixed_sessReg.F.val.F)); % replace with smallest possible number if impossible value
%     results.OLS.fixed_sessReg.F.df.pFull = pFull;
%     results.OLS.fixed_sessReg.F.df.pReduced = pReduced;
%     results.OLS.fixed_sessReg.F.df.n = n;
%     
%     % Compute and output F (considering sess1 only)
%     %Define degrees of freedom
%     pFull = length(find(any(designMatrixFull(ind_sess1,:),1)));
%     pReduced = pFull-2;
%     n = size(designMatrixFull(ind_sess1,:),1);
%     %Compute F
%     results.OLS.fixed_sessReg.Fsess1.val.F = ((RSSreduced_sess1-RSSfull_sess1)./(pFull-pReduced)) ./ (RSSfull_sess1./(n-pFull));
%     results.OLS.fixed_sessReg.Fsess1.val.F(results.OLS.fixed_sessReg.Fsess1.val.F<=0) = eps(class(results.OLS.fixed_sessReg.Fsess1.val.F)); % replace with smallest possible number if impossible value
%     results.OLS.fixed_sessReg.Fsess1.df.pFull = pFull;
%     results.OLS.fixed_sessReg.Fsess1.df.pReduced = pReduced;
%     results.OLS.fixed_sessReg.Fsess1.df.n = n;
%     
%     % Compute and output F (considering sess2 only)
%     %Define degrees of freedom
%     pFull = length(find(any(designMatrixFull(ind_sess2,:),1)));
%     pReduced = pFull-2;
%     n = size(designMatrixFull(ind_sess2,:),1);
%     %Compute F
%     results.OLS.fixed_sessReg.Fsess2.val.F = ((RSSreduced_sess2-RSSfull_sess2)./(pFull-pReduced)) ./ (RSSfull_sess2./(n-pFull));
%     results.OLS.fixed_sessReg.Fsess2.val.F(results.OLS.fixed_sessReg.Fsess2.val.F<=0) = eps(class(results.OLS.fixed_sessReg.Fsess2.val.F)); % replace with smallest possible number if impossible value
%     results.OLS.fixed_sessReg.Fsess2.df.pFull = pFull;
%     results.OLS.fixed_sessReg.Fsess2.df.pReduced = pReduced;
%     results.OLS.fixed_sessReg.Fsess2.df.n = n;
%     
%     % Compute and output F (considering only cond3)
%     %Define degrees of freedom
%     pFull = length(find(any(designMatrixFull(ind_cond3,:),1)));
%     pReduced = pFull-2;
%     n = size(designMatrixFull(ind_cond3,:),1);
%     %Compute F
%     results.OLS.fixed_sessReg.Fcond3.val.F = ((RSSreduced_cond3-RSSfull_cond3)./(pFull-pReduced)) ./ (RSSfull_cond3./(n-pFull));
%     results.OLS.fixed_sessReg.Fcond3.val.F(results.OLS.fixed_sessReg.Fcond3.val.F<=0) = eps(class(results.OLS.fixed_sessReg.Fcond3.val.F)); % replace with smallest possible number if impossible value
%     results.OLS.fixed_sessReg.Fcond3.df.pFull = pFull;
%     results.OLS.fixed_sessReg.Fcond3.df.pReduced = pReduced;
%     results.OLS.fixed_sessReg.Fcond3.df.n = n;
%     
%     % Compute and output F (considering only cond3 sess1)
%     %Define degrees of freedom
%     pFull = length(find(any(designMatrixFull(ind_cond3sess1,:),1)));
%     pReduced = pFull-2;
%     n = size(designMatrixFull(ind_cond3sess1,:),1);
%     %Compute F
%     results.OLS.fixed_sessReg.Fcond3sess1.val.F = ((RSSreduced_cond3sess1-RSSfull_cond3sess1)./(pFull-pReduced)) ./ (RSSfull_cond3sess1./(n-pFull));
%     results.OLS.fixed_sessReg.Fcond3sess1.val.F(results.OLS.fixed_sessReg.Fcond3sess1.val.F<=0) = eps(class(results.OLS.fixed_sessReg.Fcond3sess1.val.F)); % replace with smallest possible number if impossible value
%     results.OLS.fixed_sessReg.Fcond3sess1.df.pFull = pFull;
%     results.OLS.fixed_sessReg.Fcond3sess1.df.pReduced = pReduced;
%     results.OLS.fixed_sessReg.Fcond3sess1.df.n = n;
%     
%     % Compute and output F (considering only cond3 sess2)
%     %Define degrees of freedom
%     pFull = length(find(any(designMatrixFull(ind_cond3sess2,:),1)));
%     pReduced = pFull-2;
%     n = size(designMatrixFull(ind_cond3sess2,:),1);
%     %Compute F
%     results.OLS.fixed_sessReg.Fcond3sess2.val.F = ((RSSreduced_cond3sess2-RSSfull_cond3sess2)./(pFull-pReduced)) ./ (RSSfull_cond3sess2./(n-pFull));
%     results.OLS.fixed_sessReg.Fcond3sess2.val.F(results.OLS.fixed_sessReg.Fcond3sess2.val.F<=0) = eps(class(results.OLS.fixed_sessReg.Fcond3sess2.val.F)); % replace with smallest possible number if impossible value
%     results.OLS.fixed_sessReg.Fcond3sess2.df.pFull = pFull;
%     results.OLS.fixed_sessReg.Fcond3sess2.df.pReduced = pReduced;
%     results.OLS.fixed_sessReg.Fcond3sess2.df.n = n;
%     
%     clear temp temp2 sumsq good X fullModel reducedModel RSSfull RSSreduced
%     
%     %% Reformat outputs
%     results.OLS.fixed_sessReg.parameters =   reshape(results.OLS.fixed_sessReg.parameters',  [xyzsize size(results.OLS.fixed_sessReg.parameters,1)]);
%     
% %     %% Convert to percent BOLD
% %     results.OLS.fixed_sessReg.constant.brain = mean(results.OLS.fixed_sessReg.parameters(:,:,:,any(results.OLS.fixed_sessReg.designmatrixPieces.constant,1)).*(vectorlength(ones(numtime,1))/numtime),4); %vectorlength(ones(numtime,1))/numtime is the value in the design matrix
% %     con = 1./abs(results.OLS.fixed_sessReg.constant.brain) * 100;
% %     results.OLS.fixed_sessReg.parameters = bsxfun(@times,results.OLS.fixed_sessReg.parameters,con);
%     
%     %% Compute amp and delay
%     [results.OLS.fixed_sessReg.delay, results.OLS.fixed_sessReg.amp] = cart2pol(results.OLS.fixed_sessReg.parameters(:,:,:,1:2:6),results.OLS.fixed_sessReg.parameters(:,:,:,2:2:6));
%     toc
% else
%     sessModel = opt.sessModel;
% end
% 
% 
% 
% 
% 
% 
% 
% 
% tic
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MIXED-EFFECT FIT WITH SESSION EFFECT FROM FIXED-EFFECT MODEL REMOVED FROM DATA
% fprintf('*** MIXED-EFFECT FIT WITH SESSION EFFECT FROM FIXED-EFFECT MODEL REMOVED FROM DATA ***\n');
% %% Construct design matrix
% % Polynomials and motion
% poly = cell(1,numruns);
% polyNmotion = cell(1,numruns);
% constant = cell(1,numruns);
% 
% numtime = size(data{p},dimtime)*splitIn;
% for p=1:numruns/splitIn
%     tmp = constructpolynomialmatrix(numtime,0:opt.maxpolydeg(p));
%     poly((p-1)*splitIn+1:p*splitIn) = mat2cell(tmp,repmat(size(data{p},dimtime),1,splitIn),2)';
% end
% 
% numtime = size(data{p},dimtime)*splitIn;
% for p=1:numruns
%     polyNmotion{p} = cat(2,poly{p}, ...
%         opt.extraregressors{p}, ...
%         []);
%     tmpConstant = poly{p}; tmpConstant(:,2:end) = zeros(size(tmpConstant(:,2:end)));
%     constant{p} = cat(2,tmpConstant, ...
%         zeros(size(opt.extraregressors{p})), ...
%         []);
% end
% % Reformat to remerge motion regressors that were split before
% polyNmotion_b = reshape(polyNmotion,[opt.splitedIn length(polyNmotion)/opt.splitedIn])';
% constant_b = reshape(constant,[opt.splitedIn length(constant)/opt.splitedIn])';
% polyNmotion = cell(1,size(polyNmotion_b,1));
% constant = cell(1,size(polyNmotion_b,1));
% for ii = 1:size(polyNmotion_b,1)
%     polyNmotion{ii} = catcell(1,polyNmotion_b(ii,:));
%     constant{ii} = catcell(1,constant_b(ii,:));
% end
% clear polyNmotion_b constant_b
% 
% % Conditions and Sessions
% condLabel = cat(1,ones(numruns/3,1)*1,ones(numruns/3,1)*2,ones(numruns/3,1)*3);
% 
% condDesign = cell(1,3);
% polyNmotionDesign = cell(1,3);
% constantDesign = cell(1,3);
% for cond = 1:3
%     condInd = condLabel==cond;
%     condInd_polyNmotion = condInd(1:opt.splitedIn:end);
%     condDesign{cond} = blkdiag(cache.rawdesign{condInd});
%     
%     polyNmotionDesign{cond} = blkdiag(polyNmotion{condInd_polyNmotion});
%     constantDesign{cond} = blkdiag(constant{condInd_polyNmotion});
% end
% 
% % Merge conditions, sessions, polynomials and motion
% condDesign = blkdiag(condDesign{:});
% polyNmotionDesign = blkdiag(polyNmotionDesign{:});
% constantDesign = blkdiag(constantDesign{:});
% 
% results.OLS.mixed_sessRm.designmatrix = cat(2,condDesign,polyNmotionDesign);% time x regressors
% 
% % Extract some usefull params
% results.OLS.mixed_sessRm.designmatrixPieces.cond = cat(2,condDesign,zeros(size(polyNmotionDesign)));
% results.OLS.mixed_sessRm.designmatrixPieces.motion = cat(2,zeros(size(condDesign)),polyNmotionDesign);
% results.OLS.mixed_sessRm.designmatrixPieces.constant = cat(2,zeros(size(condDesign)),constantDesign);
% % close all
% % figure('WindowStyle','docked'); colormap gray
% % imagesc(results.OLS.mixed_sessRm.designmatrix)
% % figure('WindowStyle','docked'); colormap gray
% % imagesc(results.OLS.mixed_sessRm.designmatrixPieces.cond)
% % figure('WindowStyle','docked'); colormap gray
% % imagesc(results.OLS.mixed_sessRm.designmatrixPieces.motion)
% % figure('WindowStyle','docked'); colormap gray
% % imagesc(results.OLS.mixed_sessRm.designmatrixPieces.constant,[-1 1])
% 
% %% Estimate parameters with OLS
% display('computing OLS')
% results.OLS.mixed_sessRm.parameters = ...
%     mtimescell(olsmatrix2(results.OLS.mixed_sessRm.designmatrix), ...
%     cellfun(@(x) squish(x,dimdata)',data,'UniformOutput',0));  % regressors x voxels
% 
% %% Compute F for cond 3
% display('computing F stats for cond 3')
% % Set-up full and reduced design matrices
% designMatrixFull = results.OLS.mixed_sessRm.designmatrix;
% designMatrixReduced = designMatrixFull;
% designMatrixReduced(:,numruns*2/3*2+1:numruns*2) = zeros(size(designMatrixReduced(:,numruns*2/3*2+1:numruns*2)));
% 
% % Compute predictions
% fullModel = zeros(numvoxels,size(designMatrixFull,1));
% reducedModel = zeros(numvoxels,size(designMatrixReduced,1));
% for vox = 1:numvoxels
%     tmpModel = bsxfun(@times,results.OLS.mixed_sessRm.parameters(:,vox)',designMatrixFull);
%     fullModel(vox,:) = sum(tmpModel,2);
%     reducedModel(vox,:) = sum(tmpModel(:,any(designMatrixReduced,1)),2);
% end
% fullModel = reshape(fullModel,  [xyzsize size(fullModel,2)]);
% reducedModel = reshape(reducedModel,  [xyzsize size(reducedModel,2)]);
% 
% % Find indices of data to include for different F tests
% ind_cond3 = size(fullModel,4)*2/3+1:size(fullModel,4);
% ind_sess1 = [size(fullModel,4)*0/6+1:size(fullModel,4)*1/6 size(fullModel,4)*2/6+1:size(fullModel,4)*3/6 size(fullModel,4)*4/6+1:size(fullModel,4)*5/6];
% ind_sess2 = [size(fullModel,4)*1/6+1:size(fullModel,4)*2/6 size(fullModel,4)*3/6+1:size(fullModel,4)*4/6 size(fullModel,4)*5/6+1:size(fullModel,4)*6/6];
% ind_cond3sess1 = size(fullModel,4)*4/6+1:size(fullModel,4)*5/6;
% ind_cond3sess2 = size(fullModel,4)*5/6+1:size(fullModel,4)*6/6;
% 
% % Compute RSSs
% RSSfull = bsxfun(@minus, fullModel, catcell(4,data));
% RSSfull = bsxfun(@times, RSSfull, RSSfull);
% RSSfull_cond3 = sum(RSSfull(:,:,:,ind_cond3),4);
% RSSfull_sess1 = sum(RSSfull(:,:,:,ind_sess1),4);
% RSSfull_sess2 = sum(RSSfull(:,:,:,ind_sess2),4);
% RSSfull_cond3sess1 = sum(RSSfull(:,:,:,ind_cond3sess1),4);
% RSSfull_cond3sess2 = sum(RSSfull(:,:,:,ind_cond3sess2),4);
% RSSfull = sum(RSSfull,4);
% 
% RSSreduced = bsxfun(@minus, reducedModel, catcell(4,data));
% RSSreduced = bsxfun(@times, RSSreduced, RSSreduced);
% RSSreduced_cond3 = sum(RSSreduced(:,:,:,ind_cond3),4);
% RSSreduced_sess1 = sum(RSSreduced(:,:,:,ind_sess1),4);
% RSSreduced_sess2 = sum(RSSreduced(:,:,:,ind_sess2),4);
% RSSreduced_cond3sess1 = sum(RSSreduced(:,:,:,ind_cond3sess1),4);
% RSSreduced_cond3sess2 = sum(RSSreduced(:,:,:,ind_cond3sess2),4);
% RSSreduced = sum(RSSreduced,4);
% 
% % Compute and output F (considering all data)
% %Define degrees of freedom
% pFull = length(find(any(designMatrixFull,1)));
% pReduced = pFull-2;
% n = size(designMatrixFull,1);
% %Compute F
% results.OLS.mixed_sessRm.F.val.F = ((RSSreduced-RSSfull)./(pFull-pReduced)) ./ (RSSfull./(n-pFull));
% results.OLS.mixed_sessRm.F.val.F(results.OLS.mixed_sessRm.F.val.F<=0) = eps(class(results.OLS.mixed_sessRm.F.val.F)); % replace with smallest possible number if impossible value
% results.OLS.mixed_sessRm.F.df.pFull = pFull;
% results.OLS.mixed_sessRm.F.df.pReduced = pReduced;
% results.OLS.mixed_sessRm.F.df.n = n;
% 
% % Compute and output F (considering sess1 only)
% %Define degrees of freedom
% pFull = length(find(any(designMatrixFull(ind_sess1,:),1)));
% pReduced = pFull-2;
% n = size(designMatrixFull(ind_sess1,:),1);
% %Compute F
% results.OLS.mixed_sessRm.Fsess1.val.F = ((RSSreduced_sess1-RSSfull_sess1)./(pFull-pReduced)) ./ (RSSfull_sess1./(n-pFull));
% results.OLS.mixed_sessRm.Fsess1.val.F(results.OLS.mixed_sessRm.Fsess1.val.F<=0) = eps(class(results.OLS.mixed_sessRm.Fsess1.val.F)); % replace with smallest possible number if impossible value
% results.OLS.mixed_sessRm.Fsess1.df.pFull = pFull;
% results.OLS.mixed_sessRm.Fsess1.df.pReduced = pReduced;
% results.OLS.mixed_sessRm.Fsess1.df.n = n;
% 
% % Compute and output F (considering sess2 only)
% %Define degrees of freedom
% pFull = length(find(any(designMatrixFull(ind_sess2,:),1)));
% pReduced = pFull-2;
% n = size(designMatrixFull(ind_sess2,:),1);
% %Compute F
% results.OLS.mixed_sessRm.Fsess2.val.F = ((RSSreduced_sess2-RSSfull_sess2)./(pFull-pReduced)) ./ (RSSfull_sess2./(n-pFull));
% results.OLS.mixed_sessRm.Fsess2.val.F(results.OLS.mixed_sessRm.Fsess2.val.F<=0) = eps(class(results.OLS.mixed_sessRm.Fsess2.val.F)); % replace with smallest possible number if impossible value
% results.OLS.mixed_sessRm.Fsess2.df.pFull = pFull;
% results.OLS.mixed_sessRm.Fsess2.df.pReduced = pReduced;
% results.OLS.mixed_sessRm.Fsess2.df.n = n;
% 
% % Compute and output F (considering only cond3)
% %Define degrees of freedom
% pFull = length(find(any(designMatrixFull(ind_cond3,:),1)));
% pReduced = pFull-2;
% n = size(designMatrixFull(ind_cond3,:),1);
% %Compute F
% results.OLS.mixed_sessRm.Fcond3.val.F = ((RSSreduced_cond3-RSSfull_cond3)./(pFull-pReduced)) ./ (RSSfull_cond3./(n-pFull));
% results.OLS.mixed_sessRm.Fcond3.val.F(results.OLS.mixed_sessRm.Fcond3.val.F<=0) = eps(class(results.OLS.mixed_sessRm.Fcond3.val.F)); % replace with smallest possible number if impossible value
% results.OLS.mixed_sessRm.Fcond3.df.pFull = pFull;
% results.OLS.mixed_sessRm.Fcond3.df.pReduced = pReduced;
% results.OLS.mixed_sessRm.Fcond3.df.n = n;
% 
% % Compute and output F (considering only cond3 sess1)
% %Define degrees of freedom
% pFull = length(find(any(designMatrixFull(ind_cond3sess1,:),1)));
% pReduced = pFull-2;
% n = size(designMatrixFull(ind_cond3sess1,:),1);
% %Compute F
% results.OLS.mixed_sessRm.Fcond3sess1.val.F = ((RSSreduced_cond3sess1-RSSfull_cond3sess1)./(pFull-pReduced)) ./ (RSSfull_cond3sess1./(n-pFull));
% results.OLS.mixed_sessRm.Fcond3sess1.val.F(results.OLS.mixed_sessRm.Fcond3sess1.val.F<=0) = eps(class(results.OLS.mixed_sessRm.Fcond3sess1.val.F)); % replace with smallest possible number if impossible value
% results.OLS.mixed_sessRm.Fcond3sess1.df.pFull = pFull;
% results.OLS.mixed_sessRm.Fcond3sess1.df.pReduced = pReduced;
% results.OLS.mixed_sessRm.Fcond3sess1.df.n = n;
% 
% % Compute and output F (considering only cond3 sess2)
% %Define degrees of freedom
% pFull = length(find(any(designMatrixFull(ind_cond3sess2,:),1)));
% pReduced = pFull-2;
% n = size(designMatrixFull(ind_cond3sess2,:),1);
% %Compute F
% results.OLS.mixed_sessRm.Fcond3sess2.val.F = ((RSSreduced_cond3sess2-RSSfull_cond3sess2)./(pFull-pReduced)) ./ (RSSfull_cond3sess2./(n-pFull));
% results.OLS.mixed_sessRm.Fcond3sess2.val.F(results.OLS.mixed_sessRm.Fcond3sess2.val.F<=0) = eps(class(results.OLS.mixed_sessRm.Fcond3sess2.val.F)); % replace with smallest possible number if impossible value
% results.OLS.mixed_sessRm.Fcond3sess2.df.pFull = pFull;
% results.OLS.mixed_sessRm.Fcond3sess2.df.pReduced = pReduced;
% results.OLS.mixed_sessRm.Fcond3sess2.df.n = n;
% 
% clear temp temp2 sumsq good X fullModel reducedModel RSSfull RSSreduced
% 
% %% Reformat outputs
% results.OLS.mixed_sessRm.parameters =   reshape(results.OLS.mixed_sessRm.parameters',  [xyzsize size(results.OLS.mixed_sessRm.parameters,1)]);
% 
% % %% Convert to percent BOLD
% % results.OLS.mixed_sessRm.constant.brain = mean(results.OLS.mixed_sessRm.parameters(:,:,:,any(results.OLS.mixed_sessRm.designmatrixPieces.constant,1)).*(vectorlength(ones(numtime,1))/numtime),4); %vectorlength(ones(numtime,1))/numtime is the value in the design matrix
% % con = 1./abs(results.OLS.mixed_sessRm.constant.brain) * 100;
% % results.OLS.mixed_sessRm.parameters = bsxfun(@times,results.OLS.mixed_sessRm.parameters,con);
% 
% %% Compute amp and delay
% [results.OLS.mixed_sessRm.delay, results.OLS.mixed_sessRm.amp] = cart2pol(results.OLS.mixed_sessRm.parameters(:,:,:,1:2:numruns*2),results.OLS.mixed_sessRm.parameters(:,:,:,2:2:numruns*2));
% toc


end



if isfield(opt,'sessModel')
    opt.sessModel = [];
end
sessModel = [];

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




