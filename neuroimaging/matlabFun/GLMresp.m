function [results] = GLMresp(design,data,stimdur,tr,hrfmodel,hrfknobs,opt,splitIn)

%WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING
%stimdur is internally fixed to 6 in GLMestimatemodel>fitmodel_helper at 884
%WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEAL WITH INPUTS, ETC.

% input
if strcmp(hrfmodel,'SinCos')
    for i = 1:length(design)
        design{i} = repmat(design{i},1,2);
    end
elseif strcmp(hrfmodel,'resp')
    for i = 1:length(design)
        design{i} = repmat(design{i},1,stimdur/tr*2);
    end
end
if ~exist('hrfmodel','var') || isempty(hrfmodel)
    hrfmodel = 'optimize';
end
if ~exist('hrfknobs','var') || isempty(hrfknobs)
    if isequal(hrfmodel,'fir')
        hrfknobs = 20;
    elseif isequal(hrfmodel,'SinCos')
        hrfknobs = sin(2*pi*1/(stimdur/tr*2)*(0:tr:stimdur/tr*2-tr))';
        hrfknobs = normalizemax(hrfknobs);
    elseif isequal(hrfmodel,'resp')
        hrfknobs = zeros(stimdur/tr*2,1);
        hrfknobs(1) = 1;
    else
        hrfknobs = normalizemax(getcanonicalhrf(stimdur/tr,tr)');
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

if strcmp(hrfmodel,'SinCos')
    for p = 1:length(design)
        ntime = size(design{p},1);                    % number of time points
        hrfknobs(:,2) = normalizemax(cos(2*pi*1/(6*2)*(0:tr:6*2-tr))');
        for i = 1:2
            tmp(:,i) = conv2(full(design{p}(:,i)),hrfknobs(:,i));  % convolve
        end
        tmp = tmp(1:ntime,:);             % extract desired subset
        cache.rawdesign{p} = tmp; clear tmp
    end
elseif strcmp(hrfmodel,'resp')
    for p = 1:length(design)
        ntime = size(design{p},1);                    % number of time points
        hrfknobs = zeros(stimdur/tr*2);
        hrfknobs(logical(eye(size(hrfknobs)))) = 1;
        for i = 1:size(hrfknobs,2)
            tmp(:,i) = conv2(full(design{p}(:,i)),hrfknobs(:,i));  % convolve
        end
        tmp = tmp(1:ntime,:);             % extract desired subset
        cache.rawdesign{p} = tmp; clear tmp
    end
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
condLabel = cat(1,ones(numruns/3,1)*1,ones(numruns/3,1)*2,ones(numruns/3,1)*3);
sessLabel = cell2mat(opt.sessionLabel)';

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


%remove constant
results.OLS.fixed.designmatrix(:,any(results.OLS.fixed.designmatrixPieces.constant,1)) = [];
results.OLS.fixed.designmatrixPieces.constant = zeros(size(results.OLS.fixed.designmatrixPieces.constant));

%% Estimate parameters with OLS
display('computing OLS')
results.OLS.fixed.parameters = ...
    mtimescell(olsmatrix2(results.OLS.fixed.designmatrix), ...
    cellfun(@(x) squish(x,dimdata)',data,'UniformOutput',0));  % regressors x voxels

for i = 1:3*2
    respFixed(:,:,i) = results.OLS.fixed.parameters((i-1)*12+1:i*12,:);
end
results.OLS.fixed = rmfield(results.OLS.fixed,'parameters');

%percent BOLD
respFixed = (respFixed - repmat(mean(respFixed,1),[12 1 1]))./repmat(mean(respFixed,1),[12 1 1]);
respFixed = cat(4,respFixed(:,:,1:2:end),respFixed(:,:,2:2:end));

for cond = 1:3
    results.OLS.fixed.sess1.resp(:,:,:,:,cond) =   reshape(respFixed(:,:,cond,1)',  [xyzsize size(respFixed(:,:,cond,1),1)]);
    results.OLS.fixed.sess2.resp(:,:,:,:,cond) =   reshape(respFixed(:,:,cond,2)',  [xyzsize size(respFixed(:,:,cond,2),1)]);
end
clear respFixed




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MIXED-EFFECT FIT
tic
fprintf('*** MIXED-EFFECT FIT ***\n');
%% Construct design matrix
% Polynomials and motion
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
condLabel = cat(1,ones(numruns/3,1)*1,ones(numruns/3,1)*2,ones(numruns/3,1)*3);
sessLabel = cell2mat(opt.sessionLabel)';

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
    condDesign1{cond} = blkdiag(cache.rawdesign{condInd});
    
    polyNmotionDesign1{cond} = blkdiag(polyNmotion{condInd_polyNmotion});
    constantDesign1{cond} = blkdiag(constant{condInd_polyNmotion});
    
    
    sess=2;
    condInd = condLabel==cond & sessLabel==sess;
    condInd_polyNmotion = condInd(1:opt.splitedIn:end);
    condDesign2{cond} = blkdiag(cache.rawdesign{condInd});
    
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

results.OLS.mixed.designmatrix = cat(2,condDesign,polyNmotionDesign);% time x regressors
% results.OLS.mixed.designmatrixCond3 = cat(2,condDesign(end*2/3+1:end,:),polyNmotionDesign(end*2/3+1:end,:));% time x regressors

% Extract some usefull params
results.OLS.mixed.designmatrixPieces.cond = cat(2,condDesign,zeros(size(polyNmotionDesign)));
results.OLS.mixed.designmatrixPieces.motion = cat(2,zeros(size(condDesign)),polyNmotionDesign);
results.OLS.mixed.designmatrixPieces.constant = cat(2,zeros(size(condDesign)),constantDesign);
%     close all
%     figure('WindowStyle','docked'); colormap gray
%     imagesc(results.OLS.mixed.designmatrix)
%     figure('WindowStyle','docked'); colormap gray
%     imagesc(results.OLS.mixed.designmatrixPieces.cond)
%     figure('WindowStyle','docked'); colormap gray
%     imagesc(results.OLS.mixed.designmatrixPieces.motion)
%     figure('WindowStyle','docked'); colormap gray
%     imagesc(results.OLS.mixed.designmatrixPieces.constant,[-1 1])


%remove constant
results.OLS.mixed.designmatrix(:,any(results.OLS.mixed.designmatrixPieces.constant,1)) = [];
results.OLS.mixed.designmatrixPieces.constant = zeros(size(results.OLS.mixed.designmatrixPieces.constant));

%% Estimate parameters with OLS
display('computing OLS')
results.OLS.mixed.parameters = ...
    mtimescell(olsmatrix2(results.OLS.mixed.designmatrix), ...
    cellfun(@(x) squish(x,dimdata)',data,'UniformOutput',0));  % regressors x voxels

for i = 1:length(opt.sessionLabel)
    respTmp(:,:,i) = results.OLS.mixed.parameters((i-1)*12+1:i*12,:);
end
results.OLS.mixed = rmfield(results.OLS.mixed,'parameters');

respMixed1 = respTmp(:,:,cell2mat(opt.sessionLabel)==1);
respMixed1 = cat(4,respMixed1(:,:,1:end/3),respMixed1(:,:,end/3+1:end/3*2),respMixed1(:,:,end/3*2+1:end));
respMixed2 = respTmp(:,:,cell2mat(opt.sessionLabel)==2);
respMixed2 = cat(4,respMixed2(:,:,1:end/3),respMixed2(:,:,end/3+1:end/3*2),respMixed2(:,:,end/3*2+1:end));
clear respTmp

%Percent BOLD
respMixed1 = (respMixed1 - repmat(mean(respMixed1,1),[12 1 1 1]))./repmat(mean(respMixed1,1),[12 1 1 1]);
respMixed2 = (respMixed2 - repmat(mean(respMixed2,1),[12 1 1 1]))./repmat(mean(respMixed2,1),[12 1 1 1]);

for run = 1:size(respMixed1,3)
    for cond = 1:size(respMixed1,4)
        results.OLS.mixed.sess1.resp(:,:,:,:,run,cond) =   reshape(respMixed1(:,:,run,cond)',  [xyzsize size(respMixed1(:,:,run,cond),1)]);
    end
end
for run = 1:size(respMixed2,3)
    for cond = 1:size(respMixed2,4)
        results.OLS.mixed.sess2.resp(:,:,:,:,run,cond) =   reshape(respMixed2(:,:,run,cond)',  [xyzsize size(respMixed2(:,:,run,cond),1)]);
    end
end

