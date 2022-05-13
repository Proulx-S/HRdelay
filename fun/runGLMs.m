function res = runGLMs(d,p,verbose)
if ~exist('verbose','var')
    verbose = 1;
end
% function [results] = runGLM(design,data,stimdur,tr,hrfmodel,hrfknobs,opt,splitIn)


p.timeSz = zeros(size(d.data,1),1);
p.xyzSz = size(d.data{1},1:3);
p.runSz = size(d.data,1);
p.polyDeg = zeros(size(d.data,1),1);
d.poly = cell(size(d.data));
d.polyInfo = cell(size(d.data));
d.censorPts = cell(size(d.data));
for runInd = 1:size(d.data,1)
    timeSz = size(d.data{runInd},4);
    polyDeg = round(timeSz*p.tr/60/2);
    d.poly{runInd,1} = constructpolynomialmatrix(timeSz,0:polyDeg);
    d.poly{runInd,1}(:,2:end) = d.poly{runInd,1}(:,2:end) - min(d.poly{runInd,1}(:,2:end),[],1) - 1;
    d.poly{runInd,1}(:,2:end) = d.poly{runInd,1}(:,2:end) ./ max(d.poly{runInd,1}(:,2:end),[],1);
    d.polyInfo{runInd,1} = cellstr(num2str((0:polyDeg)','poly%d'))';
    p.timeSz(runInd,1) = timeSz;
    p.polyDeg(runInd,1) = polyDeg;
    d.censorPts{runInd,1} = false(timeSz,1);
    d.censorPts{runInd,1}(1:24) = true;
end
p.censorPts = d.censorPts;


%% Vein map
if verbose
    disp('Detrend, map veins and convert to %BOLD')
end
opt.excl = 0;
opt.extraOutput = 1;
% Detrend
[~,d.base,d.dataDtrd] = fitSinMixed(d,p,opt); clear opt
v.map = cell(size(d.data));
v.var = cell(size(d.data));
v.base = cell(size(d.data));
v.excl = d.excl;
for runInd = 1:size(d.data,1)
    % Vein map
    v.var{runInd} = std(d.dataDtrd{runInd}(:,:,:,~d.censorPts{runInd}),[],4);
    v.base{runInd} = d.base{runInd};
    v.map{runInd} = v.var{runInd}./v.base{runInd};
    
end

%% Convert to percent BOLD
for runInd = 1:size(d.data,1)
    d.data{runInd} = (d.data{runInd}-d.base{runInd}) ./ d.base{runInd} .* 100;
    d.dataDtrd{runInd} = d.dataDtrd{runInd} ./ d.base{runInd} .* 100;
end
    
%% 
if verbose; disp('Extract hr'); end
res1 = fitHrMixed(d,p);
if verbose; disp('Extract sin responses and fixed-effect stats'); end
res2 = fitSinMixed(d,p);

res.sin = res2.hr; res2.hr = [];
res.sinBase = res2.base; res2.base = [];
res.hr = res1.hr; res1.hr = [];
res.hrBase = res1.base; res1.base = [];
res.info = 'x X y X z X rep X cond X t';
res.sinDesign = res2.design;
res.hrDesign = res1.design; clear res1
res.infoDesign = 't X regressor';
res.featSel.F = res2.F; clear res2
res.featSel.vein.map = cat(5,...
    cat(4,v.map{d.condLabel==1 & ~d.excl}),...
    cat(4,v.map{d.condLabel==2 & ~d.excl}),...
    cat(4,v.map{d.condLabel==3 & ~d.excl})); clear v
res.dataDtrd = d.dataDtrd; clear d

function res = fitSinFixed(d,p)
%% First exclude
excl = d.excl;
if any(excl)
    rep = unique(d.repLabel(excl));
    if length(rep)>1
        error('cannot deal with multiple exclusions')
    end
    fieldList = fields(d);
    for fieldInd = 1:length(fieldList)
        d.(fieldList{fieldInd})(excl,:,:) = [];
    end
    d.repLabel(d.repLabel>rep) = d.repLabel(d.repLabel>rep)-1;
    d.runInd = (1:size(d.repLabel,1))';
    fieldList = fields(p);
    for fieldInd = 1:length(fieldList)
        if size(p.(fieldList{fieldInd}),1)==p.runSz
            p.(fieldList{fieldInd})(excl,:,:) = [];
        end
    end
    p.runSz = p.runSz - nnz(excl);
else
    d.runInd = (1:size(d.repLabel,1))';
end

%% Prepare peices of design matrix
for runInd = 1:size(d.data,1)
    t = (1:p.tr:p.tr*size(d.design{runInd},1))-1;
    d.design{runInd}(:,1) = normalizemax(sin(2*pi*1/(p.stimDur*2)*t)');
    d.design{runInd}(:,2) = normalizemax(cos(2*pi*1/(p.stimDur*2)*t)');
end
p.designInfo1 = {'sin' 'cos'};
p.designInfo2 = cellstr(num2str(sort(unique(d.condLabel)),'cond%d'))';

%% Fixed-effect
f = fitFixed(d,p);
% fieldList = fields(f);
% for fieldInd = 1:length(fieldList)
%     if isstruct(f.(fieldList{fieldInd}))
%         figure('WindowStyle','docked');
%         imagesc(catcell(1,f.(fieldList{fieldInd}).design),[-1 1])
%         title(f.(fieldList{fieldInd}).info)
%     end
% end

%% F stats
f.full = getYhat(f.full,p);
f.full = getSS(f.full,'yHat');
f.full = getYerr(f.full,d);
f.full = getSS(f.full,'yErr');


testLabel = 'act';
fullLabel = 'full';
nullLabel = 'actNull';
condInd = [1 2 3];
runInd = ismember(d.condLabel,condInd);

f.(nullLabel) = getYhat(f.(nullLabel),p);
f.(nullLabel) = getSS(f.(nullLabel),'yHat');
f.(nullLabel) = getYerr(f.(nullLabel),d);
f.(nullLabel) = getSS(f.(nullLabel),'yErr');

F.(testLabel) = getF(f.(fullLabel),f.(nullLabel),runInd);
f.(nullLabel) = rmfield(f.(nullLabel),{'yHat' 'yHatSS' 'yHatSS_n' 'yErr' 'yErrSS' 'yErrSS_n'});

res.F = F;

function [res,baseData,dataDtrd] = fitSinMixed(d,p,opt)
if ~exist('opt','var')
    opt.excl = 1;
end
if ~isfield(opt,'excl')
    opt.excl = 1;
end
if ~isfield(opt,'extraOutput')
    opt.extraOutput = 0;
end

%% First exclude
if opt.excl
    excl = d.excl;
    if any(excl)
        rep = unique(d.repLabel(excl));
        if length(rep)>1
            error('cannot deal with multiple exclusions')
        end
        fieldList = fields(d);
        for fieldInd = 1:length(fieldList)
            d.(fieldList{fieldInd})(excl,:,:) = [];
        end
        d.repLabel(d.repLabel>rep) = d.repLabel(d.repLabel>rep)-1;
        d.runInd = (1:size(d.repLabel,1))';
        fieldList = fields(p);
        for fieldInd = 1:length(fieldList)
            if size(p.(fieldList{fieldInd}),1)==p.runSz
                p.(fieldList{fieldInd})(excl,:,:) = [];
            end
        end
        p.runSz = p.runSz - nnz(excl);
    else
        d.runInd = (1:size(d.repLabel,1))';
    end
end

%% Prepare peices of design matrix
for runInd = 1:size(d.data,1)
    t = (1:p.tr:p.tr*size(d.design{runInd},1))-1;
    d.design{runInd}(:,1) = normalizemax(sin(2*pi*1/(p.stimDur*2)*t)');
    d.design{runInd}(:,2) = normalizemax(cos(2*pi*1/(p.stimDur*2)*t)');
end
p.designInfo1 = {'sin' 'cos'};
p.designInfo2 = cellstr(num2str(sort(unique(d.condLabel)),'cond%d'))';


%% Mixed-effect fit
f = fitMixed(d,p,opt);


%% F stats
if ~opt.extraOutput
    f.reduced = getRedModel(f.full);
    f.full = getYhat(f.full,p);
    f.full = getSS(f.full,'yHat');
    f.full = getYerr(f.full,d);
    f.full = getSS(f.full,'yErr');
    
    
    testLabel = 'act';
    fullLabel = 'full';
    nullLabel = 'reduced';
    condInd = [1 2 3];
    runInd = ismember(d.condLabel,condInd);
    
    f.(nullLabel) = getYhat(f.(nullLabel),p);
    f.(nullLabel) = getSS(f.(nullLabel),'yHat');
    f.(nullLabel) = getYerr(f.(nullLabel),d);
    f.(nullLabel) = getSS(f.(nullLabel),'yErr');
    
    F.(testLabel) = getF(f.(fullLabel),f.(nullLabel),runInd,1);
    f.(nullLabel) = rmfield(f.(nullLabel),{'yHat' 'yHatSS' 'yHatSS_n' 'yErr' 'yErrSS' 'yErrSS_n'});
    % imagesc(F.(testLabel).F(:,:,10))
end

%% Extra output
if opt.extraOutput
    modelName = 'full';
    fTmp = f.(modelName);
    %remove regressors of interest
    for runInd = 1:size(fTmp.design,1)
        ind = ~cellfun('isempty',fTmp.designInfo{runInd}(1,:));
        fTmp.designInfo{runInd}(:,ind) = [];
        fTmp.design{runInd}(:,ind) = [];
        fTmp.betas{runInd}(:,:,:,ind) = [];
    end
    %get detrended data from the residual of the reduced full model
    fTmp = getYhat(fTmp,p);
    fTmp = getYerr(fTmp,d);
    fTmp = rmfield(fTmp,{'yHat'});
    dataDtrd = fTmp.yErr; clear fTmp
else
    dataDtrd = [];
end

%% Extrac responses
hr = cell(size(d.data));
for runInd = 1:size(f.full.betas,1)
    designInfo = f.full.designInfo{runInd};
    hrInfo = designInfo(1,:);
    tmpSin = f.full.betas{runInd}(:,:,:,ismember(hrInfo,'sin'));
    tmpCos = f.full.betas{runInd}(:,:,:,ismember(hrInfo,'cos'));
    hr{runInd} = complex(tmpSin,tmpCos);
end
res.hr = cat(5,...
    cat(4,hr{d.condLabel==1}),...
    cat(4,hr{d.condLabel==2}),...
    cat(4,hr{d.condLabel==3}));

%% Extract baseline
baseData = cell(size(d.data));
for runInd = 1:size(f.full.betas,1)
    designInfo = f.full.designInfo{runInd};
    baseInfo = designInfo(4,:);
    baseLabel = cellstr(num2str((0:p.polyDeg(runInd))','poly%d'))';
    baseInd = ismember(baseInfo,baseLabel);
    betas = permute(f.full.betas{runInd}(:,:,:,baseInd),[4 1 2 3]);
    baseData{runInd} = nan(p.xyzSz);
    for voxInd = 1:prod(p.xyzSz)
        baseData{runInd}(voxInd) = mean(d.poly{runInd}(~d.censorPts{runInd},:)*betas(:,voxInd),1);
    end
end
res.base = cat(5,...
    cat(4,baseData{d.condLabel==1}),...
    cat(4,baseData{d.condLabel==2}),...
    cat(4,baseData{d.condLabel==3}));

res.info = 'x X y X x X rep X cond';
res.design = f.full.design{1};
res.design(f.full.censorPts{1},:) = 0;
if ~opt.extraOutput
    res.F = F;
end


function fModel = getRedModel(fModel)
for runInd = 1:size(fModel.betas,1)
    redRegrInd = cellfun('isempty',fModel.designInfo{runInd}(1,:));
    fModel.betas{runInd}(:,:,:,redRegrInd) = [];
    fModel.design{runInd}(:,redRegrInd) = [];
    fModel.designInfo{runInd}(:,redRegrInd) = [];
end
fModel.info='reduced model';





function res = fitHrMixed(d,p)
%% First exclude
excl = d.excl;
if any(excl)
    rep = unique(d.repLabel(excl));
    if length(rep)>1
        error('cannot deal with multiple exclusions')
    end
    fieldList = fields(d);
    for fieldInd = 1:length(fieldList)
        d.(fieldList{fieldInd})(excl,:,:) = [];
    end
    d.repLabel(d.repLabel>rep) = d.repLabel(d.repLabel>rep)-1;
    d.runInd = (1:size(d.repLabel,1))';
    fieldList = fields(p);
    for fieldInd = 1:length(fieldList)
        if size(p.(fieldList{fieldInd}),1)==p.runSz
            p.(fieldList{fieldInd})(excl,:,:) = [];
        end
    end
    p.runSz = p.runSz - nnz(excl);
else
    d.runInd = (1:size(d.repLabel,1))';
end

%% Prepare peices of design matrix
for runInd = 1:size(d.data,1)
    hrfknobs = zeros(p.stimDur/p.tr*2);
    hrfknobs(logical(eye(size(hrfknobs)))) = 1;
    hrfknobs(:,1) = [];
    d.design{runInd} = repmat(d.design{runInd},1,size(hrfknobs,2));
    tmp = nan(p.timeSz(runInd)+p.stimDur/p.tr*2-1,size(hrfknobs,2));
    for knobInd = 1:size(hrfknobs,2)
        tmp(:,knobInd) = conv2(full(d.design{runInd}(:,knobInd)),hrfknobs(:,knobInd));  % convolve
    end
    d.design{runInd} = tmp(1:p.timeSz(runInd,1),:); clear tmp
end
p.designInfo1 = cellstr(num2str((1:size(hrfknobs,2))','t%d'))';
p.designInfo2 = cellstr(num2str(sort(unique(d.condLabel)),'cond%d'))';

%% Mixed-effect
f = fitMixed(d,p);

%% Extract responses
hr = cell(size(d.data));
for runInd = 1:size(f.full.betas,1)
    designInfo = f.full.designInfo{runInd};
    hrInfo = designInfo(1,:);
    hr{runInd} = f.full.betas{runInd}(:,:,:,~cellfun('isempty',hrInfo));
end
res.hr = cat(6,...
    cat(5,hr{d.condLabel==1}),...
    cat(5,hr{d.condLabel==2}),...
    cat(5,hr{d.condLabel==3}));
res.hr = permute(res.hr,[1 2 3 5 6 4]);
res.hr = cat(6,zeros(size(res.hr,1:5)),res.hr);

%% Extract baseline
baseData = cell(size(d.data));
for runInd = 1:size(f.full.betas,1)
    designInfo = f.full.designInfo{runInd};
    baseInfo = designInfo(4,:);
    baseLabel = cellstr(num2str((0:p.polyDeg(runInd))','poly%d'))';
    baseInd = all(ismember(baseInfo,baseLabel),1);
    baseData{runInd} = mean(f.full.betas{runInd}(:,:,:,baseInd),4);
end
baseData = cat(5,...
    cat(4,baseData{d.condLabel==1}),...
    cat(4,baseData{d.condLabel==2}),...
    cat(4,baseData{d.condLabel==3}));

res.hr = res.hr + baseData;
res.base = baseData;
res.info = 'x X y X x X rep X cond X t';
res.design = f.full.design{1};
res.design(f.full.censorPts{1},:) = 0;

% tmp = permute(res.hr,[6 4 5 1 2 3]);
% tmp = mean(tmp(:,:,:,:),4);
% plot(tmp(:,:))


function res = runFit(d,p,opt)
%% First exclude
excl = d.excl;
rep = unique(d.repLabel(excl));
if length(rep)>1
    error('cannot deal with multiple exclusions')
end
fieldList = fields(d);
for fieldInd = 1:length(fieldList)
    d.(fieldList{fieldInd})(excl,:,:) = [];
end
d.repLabel(d.repLabel>rep) = d.repLabel(d.repLabel>rep)-1;
d.runInd = (1:size(d.repLabel,1))';
fieldList = fields(p);
for fieldInd = 1:length(fieldList)
    if size(p.(fieldList{fieldInd}),1)==p.runSz
        p.(fieldList{fieldInd})(excl,:,:) = [];
    end
end
p.runSz = p.runSz - nnz(excl);

%% Prepare peices of design matrix
switch opt.hrf
    case 'sin'
        for runInd = 1:size(d.data,1)
            t = (1:p.tr:p.tr*size(d.design{runInd},1))-1;
            d.design{runInd}(:,1) = normalizemax(sin(2*pi*1/(p.stimDur*2)*t)');
            d.design{runInd}(:,2) = normalizemax(cos(2*pi*1/(p.stimDur*2)*t)');
        end
        p.designInfo1 = {'sin' 'cos'};
        p.designInfo2 = cellstr(num2str(sort(unique(d.condLabel)),'cond%d'))';
    case 'hr'
        for runInd = 1:size(d.data,1)
            hrfknobs = zeros(p.stimDur/p.tr*2);
            hrfknobs(logical(eye(size(hrfknobs)))) = 1;
            hrfknobs(:,1) = [];
            d.design{runInd} = repmat(d.design{runInd},1,size(hrfknobs,2));
            tmp = nan(p.timeSz(runInd)+p.stimDur/p.tr*2-1,size(hrfknobs,2));
            for knobInd = 1:size(hrfknobs,2)
                tmp(:,knobInd) = conv2(full(d.design{runInd}(:,knobInd)),hrfknobs(:,knobInd));  % convolve
            end
            d.design{runInd} = tmp(1:p.timeSz(runInd,1),:); clear tmp
        end
        p.designInfo1 = cellstr(num2str((1:size(hrfknobs,2))','t%d'))';
        p.designInfo2 = cellstr(num2str(sort(unique(d.condLabel)),'cond%d'))';
    otherwise
        error('X')
end

%% Mixed-effect
disp('Mixed-Effect')
f = fitMixed(d,p,opt);
info = squeeze(f.full.designInfo(:,:,:,:,:,1))';
hrInd = all(~cellfun('isempty',info),1);
info = squeeze(f.full.designInfo(:,:,:,:,:,4))';
baseInd = ismember(info(1,:),'poly0');

base = nan([p.xyzSz size(f.full.betas,5)]);
for runInd = 1:size(f.full.betas,5)
    base(:,:,:,runInd) = f.full.betas(:,:,:,baseInd,runInd);
end
hr = repmat(base,[1 1 1 1 p.stimDur*2./p.tr]);
for runInd = 1:size(f.full.betas,5)
    hr(:,:,:,runInd,2:end) = hr(:,:,:,runInd,2:end) + permute(f.full.betas(:,:,:,hrInd,runInd),[1 2 3 5 4]);
end


%% Fixed-effect
disp('Fixed-Effect')
f = fitFixed(d,p,opt);
fieldList = fields(f);
for fieldInd = 1:length(fieldList)
    if isstruct(f.(fieldList{fieldInd}))
        figure('WindowStyle','docked');
        imagesc(catcell(1,f.(fieldList{fieldInd}).design))
        title(f.(fieldList{fieldInd}).info)
    end
end

%% Extract resp and brain
[resp,~] = getBetas(f,p);
tmp = permute(resp.hr,[4 5 1 2 3]);
tmp = mean(tmp(:,:,:),3)';
plot([0 0 0; tmp])



%% F stats
f.full = getYhat(f.full,p);
f.full = getSS(f.full,'yHat');
f.full = getYerr(f.full,d);
f.full = getSS(f.full,'yErr');


testLabel = 'act';
fullLabel = 'full';
nullLabel = 'actNull';
condInd = [1 2 3];
runInd = ismember(d.condLabel,condInd);

f.(nullLabel) = getYhat(f.(nullLabel),p);
f.(nullLabel) = getSS(f.(nullLabel),'yHat');
f.(nullLabel) = getYerr(f.(nullLabel),d);
f.(nullLabel) = getSS(f.(nullLabel),'yErr');

F.(testLabel) = getF(f.(fullLabel),f.(nullLabel),runInd);
f.(nullLabel) = rmfield(f.(nullLabel),{'yHat' 'yHatSS' 'yHatSS_n'});


testLabel = 'cond1v2v3';
fullLabel = 'full';
nullLabel = 'cond1v2v3null';
condInd = [1 2 3];
runInd = ismember(d.condLabel,condInd);

f.(nullLabel) = getYhat(f.(nullLabel),p);
f.(nullLabel) = getSS(f.(nullLabel),'yHat');
f.(nullLabel) = getYerr(f.(nullLabel),d);
f.(nullLabel) = getSS(f.(nullLabel),'yErr');

F.(testLabel) = getF(f.(fullLabel),f.(nullLabel),runInd);
f.(nullLabel) = rmfield(f.(nullLabel),{'yHat' 'yHatSS' 'yHatSS_n'});


testLabel = 'cond1v2';
fullLabel = 'full';
nullLabel = 'cond1v2null';
condInd = [1 2];
runInd = ismember(d.condLabel,condInd);

f.(nullLabel) = getYhat(f.(nullLabel),p);
f.(nullLabel) = getSS(f.(nullLabel),'yHat');
f.(nullLabel) = getYerr(f.(nullLabel),d);
f.(nullLabel) = getSS(f.(nullLabel),'yErr');

F.(testLabel) = getF(f.(fullLabel),f.(nullLabel),runInd);
f.(nullLabel) = rmfield(f.(nullLabel),{'yHat' 'yHatSS' 'yHatSS_n'});


testLabel = 'cond1v3';
fullLabel = 'full';
nullLabel = 'cond1v3null';
condInd = [1 3];
runInd = ismember(d.condLabel,condInd);

f.(nullLabel) = getYhat(f.(nullLabel),p);
f.(nullLabel) = getSS(f.(nullLabel),'yHat');
f.(nullLabel) = getYerr(f.(nullLabel),d);
f.(nullLabel) = getSS(f.(nullLabel),'yErr');

F.(testLabel) = getF(f.(fullLabel),f.(nullLabel),runInd);
f.(nullLabel) = rmfield(f.(nullLabel),{'yHat' 'yHatSS' 'yHatSS_n'});


testLabel = 'cond2v3';
fullLabel = 'full';
nullLabel = 'cond2v3null';
condInd = [2 3];
runInd = ismember(d.condLabel,condInd);

f.(nullLabel) = getYhat(f.(nullLabel),p);
f.(nullLabel) = getSS(f.(nullLabel),'yHat');
f.(nullLabel) = getYerr(f.(nullLabel),d);
f.(nullLabel) = getSS(f.(nullLabel),'yErr');

F.(testLabel) = getF(f.(fullLabel),f.(nullLabel),runInd);
f.(nullLabel) = rmfield(f.(nullLabel),{'yHat' 'yHatSS' 'yHatSS_n'});


function f = getSS(f,field)
f.([field 'SS']) = cell(size(f.(field)));
f.([field 'SS_n']) = nan(size(f.(field)));
for runInd = 1:size(f.(field),1)
    f.([field 'SS']){runInd} = sum(f.(field){runInd}(:,:,:,~f.censorPts{runInd}).^2,4);
    f.([field 'SS_n'])(runInd) = nnz(~f.censorPts{runInd});
end


function f = getYerr(f,d)
f.yErr = cell(size(f.yHat));
for runInd = 1:size(f.yHat,1)
    f.yErr{runInd} = d.data{runInd} - f.yHat{runInd};
end

function f = getYhat(f,p)
f.yHat = cell(size(f.design,1),1);
if iscell(f.betas)
    for runInd = 1:size(f.design,1)
        betas = permute(f.betas{runInd},[4 1 2 3]);
        f.yHat{runInd} = nan([p.timeSz(runInd) p.xyzSz]);
        for voxInd = 1:prod(p.xyzSz)
            f.yHat{runInd}(:,voxInd) = sum(f.design{runInd}.*betas(:,voxInd)',2);
        end
        f.yHat{runInd} = permute(f.yHat{runInd},[2 3 4 1]);
    end
else
    betas = permute(f.betas,[4 1 2 3]);
    for runInd = 1:size(f.design,1)
        f.yHat{runInd} = nan([p.timeSz(runInd) p.xyzSz]);
        for voxInd = 1:prod(p.xyzSz)
            f.yHat{runInd}(:,voxInd) = sum(f.design{runInd}.*betas(:,voxInd)',2);
        end
        f.yHat{runInd} = permute(f.yHat{runInd},[2 3 4 1]);
    end
end

function F = getF(fFull,fRed,runInd,randEffectFlag)
if ~exist('randEffectFlag','var')
    randEffectFlag = false;
end
% Process design matrix just for explicit outputs
% censor points
for i = 1:length(runInd)
    fFull.design{i}(fFull.censorPts{i},:) = 0;
    fRed.design{i}(fRed.censorPts{i},:) = 0;
end
% censor runs
for i = 1:length(runInd)
    if ~runInd(i)
        fFull.design{i}(:) = 0;
        fRed.design{i}(:) = 0;    
    end
end

censorPts = catcell(1,fFull.censorPts(runInd));
if randEffectFlag
    designFull = blkdiag(fFull.design{runInd});
    designFull = designFull(~censorPts,:);
    designRed = blkdiag(fRed.design{runInd});
    designRed = designRed(~censorPts,:);
else
    designFull = catcell(1,fFull.design(runInd));
    designFull = designFull(~censorPts,:);
    designRed = catcell(1,fRed.design(runInd));
    designRed = designRed(~censorPts,:);
end

pFull = nnz(any(designFull,1));
pRed = nnz(any(designRed,1));
n = sum(fFull.yErrSS_n(runInd));
if n~=sum(fFull.yHatSS_n(runInd))
    error('X')
end

yErrSSfull = sum(catcell(4,fFull.yErrSS(runInd)),4);
yErrSSred = sum(catcell(4,fRed.yErrSS(runInd)),4);
% yHatSSfull = sum(catcell(4,fFull.yHatSS(runInd)),4);
% yHatSSred = sum(catcell(4,fRed.yHatSS(runInd)),4);

% The following two methods are equivalent, the second demands (a little) less
% computations (no need to compute the residuals of the reduced model)
F.F  = ( (yErrSSred-yErrSSfull)./ (pFull-pRed) ) ./ ( yErrSSfull./(n-pFull) );
% F.F = ( (yHatSSfull-yHatSSred)./ (pFull-pRed) ) ./ ( yErrSSfull./(n-pFull) );

F.p = fcdf(F.F,pFull-pRed,n-pFull,'upper');
F.df.pFull = pFull;
F.df.pRed = pRed;
F.df.n = n;
F.design.full = catcell(1,fFull.design);
F.design.red = catcell(1,fRed.design);



function r = getResidualSS(r,d)
r.rss = cell(size(r.r));
r.n = nan(size(r.r));
for runInd = 1:length(r.r)
    rs = r.r{runInd}.^2;
    r.rss{runInd} = sum(rs(:,:,:,~d.censorPts{runInd}),4);
    r.n(runInd) = nnz(~d.censorPts{runInd});
end



function r = getResidual(fitRes,d,p)
betas = permute(fitRes.betas,[4 1 2 3]);
data = permute(catcell(4,d.data),[4 1 2 3]);
design = catcell(1,fitRes.design);
r.r = zeros(size(data));
for voxInd = 1:prod(p.xyzSz)
    yHat = sum(design.*betas(:,voxInd)',2);
    r.r(:,voxInd) = data(:,voxInd) - yHat;
end
r.r = permute(mat2cell(permute(r.r,[2 3 4 1]),p.xyzSz(1),p.xyzSz(2),p.xyzSz(3),p.timeSz),[4 1 2 3]);

function [resp,base] = getBetas(fitRes,p,opt)
if ~exist('opt','var')
    opt.getBase = 0;
end
reg1list = unique(fitRes.full.designInfo(1,:));
reg1list = reg1list(~cellfun('isempty',reg1list));
reg2List = unique(fitRes.full.designInfo(2,:));
reg2List = reg2List(~cellfun('isempty',reg2List));
if opt.getBase
    indBase = ismember(fitRes.full.designInfo(4,:),'poly0');
end


resp = nan([p.xyzSz length(reg2List) length(reg1list)]);
if opt.getBase
    base = nan([p.xyzSz length(reg2List) 1]);
else
    base = [];
end
for reg1ind = 1:length(reg1list)
    for reg2Ind = 1:length(reg2List)
        ind = ismember(fitRes.full.designInfo(1,:),reg1list{reg1ind}) & ismember(fitRes.full.designInfo(2,:),reg2List{reg2Ind});
        resp(:,:,:,reg2Ind,reg1ind) = fitRes.full.betas(:,:,:,ind);
        if opt.getBase
            error('code that')
        end
    end
end
switch fitRes.info
    case 'hr'
        tmp = resp; clear resp
        resp.hr = tmp;
        resp.info = 'hr: X x Y x Z x cond x time';
    case 'sin'
        tmp = complex(resp(:,:,:,:,1),resp(:,:,:,:,2)); clear resp
        resp.sin = tmp;
        resp.info = 'sin: X x Y x Z x cond';
    otherwise
        error('X')
end

function [res,poly2] = fitMixed(d,p,opt)
if ~exist('opt','var')
    opt.extraOutput = 0;
end
% if ~isfield(opt,'extraOutput')
%     opt.extraOutput = 0;
% end
% design
[designFull,designFullInfo] = getDesign(d,p);
% Polynomial regressors
[poly,polyInfo] = getPoly(d,p);
% Extra regressors (e.g. motion)
if p.doMotion
    motionInfo = [repmat({''},2,length(motionInfo)); motionInfo; repmat({''},1,length(motionInfo))];
else
    motion = cell(size(designFull));
    motionInfo = {};
end

modelName = 'full';
modelLabel = 'full model';
design = designFull;
designInfo = cell(size(design));
tmpDesign = cat(2,design,motion,poly);
res.(modelName).betas = cell(size(design));
res.(modelName).design = cell(size(design));
res.(modelName).censorPts = cell(size(design));
res.(modelName).designInfo = cell(size(design));
for runInd = 1:size(design,1)
    design{runInd} = cat(2,tmpDesign{runInd,:});
    ind = any(design{runInd},1);
    design{runInd} = design{runInd}(:,ind);
    tmp = cat(2,designFullInfo,motionInfo,polyInfo);
    designInfo{runInd} = tmp(:,ind);
    tmp = computeOLS(d.data(runInd),design(runInd),designInfo{runInd},p.censorPts(runInd));
    
    res.(modelName).betas{runInd} = tmp.betas;
    res.(modelName).design(runInd) = tmp.design(1);
    res.(modelName).censorPts(runInd) = tmp.censorPts(1);
    res.(modelName).designInfo{runInd} = tmp.designInfo;
end
res.(modelName).info = modelLabel;

% % Model 2
% modelName = 'actNull';
% modelLabel = 'Activation null model';
% tmp = cat(2,motion,poly);
% for runInd = 1:size(design,1)
%     design{runInd} = cat(2,tmp{runInd,:});
% end
% clear tmp
% designInfo = cat(2,motionInfo,polyInfo);
% res.(modelName) = computeOLS(d.data,design,designInfo,p.censorPts);
% res.(modelName).info = modelLabel;
% % imagesc(catcell(1,res.(modelName).design))
% % res.(modelName).designInfo'



% if opt.extraOutput
%     modelName = 'null';
%     modelLabel = 'null model';
%     design = designFull;
%     designInfo = cell(size(design));
%     tmpDesign = cat(2,motion,poly);
%     res.(modelName).betas = cell(size(design));
%     res.(modelName).design = cell(size(design));
%     res.(modelName).censorPts = cell(size(design));
%     res.(modelName).designInfo = cell(size(design));
%     for runInd = 1:size(design,1)
%         design{runInd} = cat(2,tmpDesign{runInd,:});
%         ind = any(design{runInd},1);
%         design{runInd} = design{runInd}(:,ind);
%         tmp = cat(2,motionInfo,polyInfo);
%         designInfo{runInd} = tmp(:,ind);
%         tmp = computeOLS(d.data(runInd),design(runInd),designInfo{runInd},p.censorPts(runInd));
%         
%         res.(modelName).betas{runInd} = tmp.betas;
%         res.(modelName).design(runInd) = tmp.design(1);
%         res.(modelName).censorPts(runInd) = tmp.censorPts(1);
%         res.(modelName).designInfo{runInd} = tmp.designInfo;
%     end
%     res.(modelName).info = modelLabel;
% end


function res = fitFixed(d,p)
% design
[designFull,designFullInfo] = getDesign(d,p);
% Polynomial regressors
[poly,polyInfo] = getPoly(d,p);
% Extra regressors (e.g. motion)
if p.doMotion
    motionInfo = [repmat({''},2,length(motionInfo)); motionInfo; repmat({''},1,length(motionInfo))];
else
    motion = cell(size(d.data));
    motionInfo = {};
end

%% Model 1
modelName = 'full';
modelLabel = 'full model';
design = designFull;
designInfo = designFullInfo;
tmp = cat(2,design,motion,poly);
for runInd = 1:size(design,1)
    design{runInd} = cat(2,tmp{runInd,:});
end
clear tmp
designInfo = cat(2,designInfo,motionInfo,polyInfo);
res.(modelName) = computeOLS(d.data,design,designInfo,p.censorPts);
res.(modelName).info = modelLabel;
% imagesc(catcell(1,res.(modelName).design))
% res.(modelName).designInfo'

%% Model 2
modelName = 'actNull';
modelLabel = 'Activation null model';
tmp = cat(2,motion,poly);
for runInd = 1:size(design,1)
    design{runInd} = cat(2,tmp{runInd,:});
end
clear tmp
designInfo = cat(2,motionInfo,polyInfo);
res.(modelName) = computeOLS(d.data,design,designInfo,p.censorPts);
res.(modelName).info = modelLabel;
% imagesc(catcell(1,res.(modelName).design))
% res.(modelName).designInfo'

%% Model 3
modelName = 'cond1v2v3null';
modelLabel = 'cond1v2v3-null model';
condToMerge = {'cond1' 'cond2' 'cond3'};
condToMergeLabel = {'cond1+2+3'};
[design,designInfo] = mergeCond(designFull,designFullInfo,condToMerge,condToMergeLabel);
tmp = cat(2,design,motion,poly);
for runInd = 1:size(design,1)
    design{runInd} = cat(2,tmp{runInd,:});
end
clear tmp
designInfo = cat(2,designInfo,motionInfo,polyInfo);
res.(modelName) = computeOLS(d.data,design,designInfo,p.censorPts);
res.(modelName).info = modelLabel;
% imagesc(catcell(1,res.(modelName).design))
% res.(modelName).designInfo'

%% Model 4
modelName = 'cond1v2null';
modelLabel = 'cond1v2-null model';
condToMerge = {'cond1' 'cond2'};
condToMergeLabel = {'cond1+2'};
[design,designInfo] = mergeCond(designFull,designFullInfo,condToMerge,condToMergeLabel);
tmp = cat(2,design,motion,poly);
for runInd = 1:size(design,1)
    design{runInd} = cat(2,tmp{runInd,:});
end
clear tmp
designInfo = cat(2,designInfo,motionInfo,polyInfo);
res.(modelName) = computeOLS(d.data,design,designInfo,p.censorPts);
res.(modelName).info = modelLabel;
% imagesc(catcell(1,res.(modelName).design))
% res.(modelName).designInfo'


%% Model 5
modelName = 'cond1v3null';
modelLabel = 'cond1v3-null model';
condToMerge = {'cond1' 'cond3'};
condToMergeLabel = {'cond1+3'};
[design,designInfo] = mergeCond(designFull,designFullInfo,condToMerge,condToMergeLabel);
tmp = cat(2,design,motion,poly);
for runInd = 1:size(design,1)
    design{runInd} = cat(2,tmp{runInd,:});
end
clear tmp
designInfo = cat(2,designInfo,motionInfo,polyInfo);
res.(modelName) = computeOLS(d.data,design,designInfo,p.censorPts);
res.(modelName).info = modelLabel;
% imagesc(catcell(1,res.(modelName).design))
% res.(modelName).designInfo'

%% Model 6
modelName = 'cond2v3null';
modelLabel = 'cond2v3-null model';
condToMerge = {'cond2' 'cond3'};
condToMergeLabel = {'cond2+3'};
[design,designInfo] = mergeCond(designFull,designFullInfo,condToMerge,condToMergeLabel);
tmp = cat(2,design,motion,poly);
for runInd = 1:size(design,1)
    design{runInd} = cat(2,tmp{runInd,:});
end
clear tmp
designInfo = cat(2,designInfo,motionInfo,polyInfo);
res.(modelName) = computeOLS(d.data,design,designInfo,p.censorPts);
res.(modelName).info = modelLabel;
% imagesc(catcell(1,res.(modelName).design))
% res.(modelName).designInfo'

function [design,designInfo] = getDesign(d,p)

condList = sort(unique(d.condLabel));
tmp1 = d.design;
runLabel = cell(size(d.design));
for runInd = 1:size(tmp1,1)
    tmp1{runInd,1} = zeros(size(tmp1{runInd,1}));
    runLabel{runInd} = ones(size(tmp1{runInd,1},1),1).*runInd;
end
design = cell(1,length(condList));
designInfo = cell(1,length(condList));
for condInd = 1:length(condList)
    cond = condList(condInd);
    tmp2 = tmp1;
    tmp2(d.condLabel==cond) = d.design(d.condLabel==cond);
    design{condInd} = catcell(1,tmp2);
    designInfo{condInd} = repmat(p.designInfo2(condInd),[1 size(design{condInd},2)]);
end
clear tmp1 tmp2
design = catcell(2,design);
designInfo1 = repmat(p.designInfo1,[1 size(designInfo,2)]);
designInfo2 = catcell(2,designInfo);
designInfo = cat(1,designInfo1,designInfo2); clear designInfo1 designInfo2
designInfo = cat(1,designInfo,repmat({''},2,length(designInfo)));
design = mat2cell(design,p.timeSz,size(design,2));

function [poly,polyInfo] = getPoly(d,p)

poly = blkdiag(d.poly{:});
polyInfo = catcell(1,d.polyInfo)';
polyInfo = polyInfo(:)';
polyInfo = [repmat({''},3,length(polyInfo)); polyInfo; repmat({''},0,length(polyInfo))];
poly = mat2cell(poly,p.timeSz,size(poly,2));

function [design,designInfo] = mergeCond(design,designInfo,condToMerge,condToMergeLabel)
% Merge conditions to merge
condToMergeInd = ismember(designInfo(2,:),condToMerge);
regList = unique(designInfo(1,:));
designMerged = cell(size(design));
designNonMerged = cell(size(design));
designTmp = cell(size(design));
for runInd = 1:size(design,1)
    designMerged{runInd} = nan([size(design{runInd,1},1) size(regList,2)]);
    for regInd = 1:length(regList)
        ind = ismember(designInfo(1,:),regList{regInd});
        ind = ind & condToMergeInd;
        designMerged{runInd}(:,regInd) = sum(design{runInd}(:,ind),2);
    end
    ind = ~ismember(designInfo(2,:),condToMerge);
    designNonMerged{runInd} = design{runInd}(:,ind);
    designTmp{runInd} = cat(2,designMerged{runInd},designNonMerged{runInd});
end
design = designTmp; clear designTmp
ind = ~ismember(designInfo(2,:),condToMerge);
designInfoNonMerged = designInfo(:,ind);
designInfoMerged = repmat({''},[size(designInfo,1) size(regList,2)]);
for regInd = 1:length(regList)
    designInfoMerged(1,regInd) = regList(regInd);
    designInfoMerged(2,regInd) = condToMergeLabel;
end
designInfo = cat(2,designInfoMerged,designInfoNonMerged); clear designInfoMerged designInfoNonMerged

function res = computeOLS(X,design,designInfo,censorPts)
% censor
designUncensored = design;
for runInd = 1:length(X)
    X{runInd}(:,:,:,censorPts{runInd}) = [];
    design{runInd}(censorPts{runInd},:) = [];
end
xyzSz = size(X{1},[1 2 3]);
res.betas = ...
    mtimescell(olsmatrix(design), ...
    cellfun(@(x) squish(x,3)',X,'UniformOutput',0));  % regressors x voxels
res.betas = permute(reshape(res.betas,[size(res.betas,1) xyzSz]),[2 3 4 1]);
res.design = designUncensored;
res.censorPts = censorPts;
res.designInfo = designInfo;




function f = constructpolynomialmatrix(n,degrees)
f = [];
temp = linspace(-1,1,n)';
for p=1:length(degrees)
    % construct polynomial
    polyvector = temp .^ degrees(p);

    % orthogonalize with respect to earlier polynomials
    polyvector = projectionmatrix(f)*polyvector;

    % record
    f = cat(2,f,polyvector);
end


function f = projectionmatrix(X)

% function f = projectionmatrix(X)
%
% <X> is samples x parameters
%
% what we want to do is to perform a regression using <X>
% and subtract out the fit.  this is accomplished by
% y-X*inv(X'*X)*X'*y = (I-X*inv(X'*X)*X')*y = f*y
% where y is the data (samples x cases).
%
% what this function does is to return <f> which has
% dimensions samples x samples.  to accomplish this,
% we rely heavily on olsmatrix.m.
%
% if <X> has no parameters, the output of this function is 1.
%
% history:
% - 2013/08/18 - in the cases of empty <X>, we now return 1 instead of [].
%
% example:
% x = sort(randn(100,1));
% x2 = projectionmatrix(constructpolynomialmatrix(100,0:1))*x;
% figure; hold on; plot(x,'r-'); plot(x2,'g-');

% handle corner case
if isempty(X)
    f = 1;
    return;
end

% do it
f = eye(size(X,1)) - X*olsmatrix(X);



function f = olsmatrix(X,mode)
if iscell(X)
    X = catcell(1,X);
end

% function f = olsmatrix(X,mode)
%
% <X> is samples x parameters
% <mode> (optional) is
%   0 means normal operation
%   1 means use inv instead of \ and omit unit-length normalization.
%     the point of this mode is to reduce memory usage (i think).
%   default: 0.
%
% what we want to do is to perform OLS regression using <X>
% and obtain the parameter estimates.  this is accomplished
% by inv(X'*X)*X'*y = f*y where y is the data (samples x cases).
%
% what this function does is to return <f> which has dimensions
% parameters x samples. 
%
% to ensure well-conditioning, we unit-length normalize each column
% of <X> before doing the inv.m operation.  also, we actually use
% left-slash (\), which apparently is equivalent to inv.m but faster
% and more accurate (see code for details).  if you pass <mode> as 1,
% we omit the normalization and use the inv method instead of the \ method.
%
% also, in the case that one or more regressors in <X> are all zero, then
% without special handling, then this case will result in warnings and
% NaN results.  what we do is to explicitly ensure that all-zero regressors
% are ignored and that there are all-zero rows for these regressors 
% in <f>.  this makes it such that the weights estimated for these 
% regressors are simply zero.
%
% history:
% 2011/06/27 - explicitly ignore 0
% 2011/03/29 - mode 1 now omits unit length normalization.
% 2011/03/28 - add <mode> input.
% 2010/08/05: now, we do not try to detect weird cases with low variance.
%             instead, we just blindly unit-length normalize each column.

% input
if ~exist('mode','var') || isempty(mode)
  mode = 0;
end

% deal with degenerate regressors
good = ~all(X==0,1);

% initialize result
f = zeros(size(X,2),size(X,1));

% do it
switch mode
case 0
  [X,len] = unitlength(X(:,good),1,[],0);
  temp = diag(1./len)*((X'*X)\X');  % hm, suggested by M-Lint
case 1
  temp = inv(X(:,good)'*X(:,good))*X(:,good)';
end

% return
if any(good)
  f(good,:) = temp;
end



function [f,len,sc] = unitlength(m,dim,flag,wantcaution,sc)

% function [f,len,sc] = unitlength(m,dim,flag,wantcaution,sc)
%
% <m> is a matrix
% <dim> (optional) is the dimension of interest.
%   if supplied, normalize each case oriented along <dim> to have unit length.
%   if [] or not supplied, normalize length globally.
% <flag> (optional) is
%   0 means normal case
%   1 means make length sqrt(n) where n is number of non-NaN entries
%   default: 0
% <wantcaution> (optional) is whether to perform special handling of 
%   weird cases where the length of <m> is very small to start with (see zerodiv.m).
%   default: 1.
% <sc> (optional) is a special case.  supply this and we will use it instead of
%   calculating the actual scale factor.  also, <len> will be returned as [].
%   to indicate not-supplied, pass in [].
%
% unit-length normalize <m> via scaling, operating either on individual cases or globally.
% the output <f> has the same dimensions as <m>.  also, return <len> which is
% the vector length of <m> along <dim>.  when <dim> is [], <len> is a scalar;
% otherwise, <len> is the same dimensions as <m> except collapsed along <dim>.
% also, return <sc> which is the scale factor divided from <m>.  the dimensions
% of <sc> is the same as <len>.
%
% we ignore NaNs gracefully.
%
% note some weird cases:
%   unitlength([]) is [].
%   unitlength([0 0]) is [NaN NaN].
%   unitlength([NaN NaN]) is [NaN NaN].
%
% history:
% 2014/04/27 - oops, make sure NaN is casted to class of <m>
% 2011/06/27 - oops. handle empty case explicitly (it would have crashed)
%
% example:
% a = [3 0 NaN];
% isequalwithequalnans(unitlength(a),[1 0 NaN])

% input
if ~exist('dim','var') || isempty(dim)
  dim = [];
end
if ~exist('flag','var') || isempty(flag)
  flag = 0;
end
if ~exist('wantcaution','var') || isempty(wantcaution)
  wantcaution = 1;
end

% handle degenerate case up front
if isempty(m)
  f = [];
  len = [];
  sc = [];
  return;
end
  
% figure out len and sc
if ~exist('sc','var') || isempty(sc)

  % figure out vector length
  len = vectorlength(m,dim);
  
  % figure out scale factor
  if flag==1
    if isempty(dim)
      temp = sqrt(sum(~isnan(m(:))));
    else
      temp = sqrt(sum(~isnan(m),dim));
    end
    sc = len./temp;
  else
    sc = len;
  end

else
  len = [];
end

% ok, do it
f = bsxfun(@(x,y) zerodiv(x,y,cast(NaN,class(m)),wantcaution),m,sc);


% HM, IS THIS SLOWER OR FASTER:
% if isempty(dim)
%   f = zerodiv(m,sc,NaN,wantcaution);
% else
%   f = zerodiv(m,repmat(sc,copymatrix(ones(1,ndims(m)),dim,size(m,dim))),NaN,wantcaution);
% end


function f = vectorlength(m,dim)

% function f = vectorlength(m,dim)
%
% <m> is a matrix
% <dim> (optional) is the dimension of interest.
%   if supplied, calculate vector length of each case oriented along <dim>.
%   if [] or not supplied, calculate vector length of entire matrix
%
% calculate vector length of <m>, either of individual cases (in which case
% the output is the same as <m> except collapsed along <dim>) or globally
% (in which case the output is a scalar).
%
% we ignore NaNs gracefully.
% 
% note some weird cases:
%   vectorlength([]) is [].
%   vectorlength([NaN NaN]) is 0
%
% example:
% a = [1 1];
% isequal(vectorlength(a),sqrt(2))
% a = [1 NaN; NaN NaN];
% isequal(vectorlength(a,1),[1 0])

% deal with NaNs
m(isnan(m)) = 0;

% handle weird case up front
if isempty(m)
  f = [];
  return;
end

% do it
if ~exist('dim','var') || isempty(dim)
  f = sqrt(dot(m(:),m(:),1));
else
  f = sqrt(dot(m,m,dim));
end



function f = zerodiv(x,y,val,wantcaution)

% function f = zerodiv(x,y,val,wantcaution)
% 
% <x>,<y> are matrices of the same size or either or both can be scalars.
% <val> (optional) is the value to use when <y> is 0.  default: 0.
% <wantcaution> (optional) is whether to perform special handling of weird
%   cases (see below).  default: 1.
%
% calculate x./y but use <val> when y is 0.
% if <wantcaution>, then if the absolute value of one or more elements of y is 
%   less than 1e-5 (but not exactly 0), we issue a warning and then treat these 
%   elements as if they are exactly 0.
% if not <wantcaution>, then we do nothing special.
%
% note some weird cases:
%   if either x or y is [], we return [].
%   NaNs in x and y are handled in the usual way.
%
% history:
% 2011/02/02 - in the case that y is not a scalar and wantcaution is set to 0, 
%              we were allowing division by 0 to result in Inf and *not* replaced 
%              with val as desired.  big mistake.  we have now fixed this.
%
% example:
% isequalwithequalnans(zerodiv([1 2 3],[1 0 NaN]),[1 0 NaN])

% input
if nargin < 4  % need optimal speed so try to bypass in the fully specified case if we can
  if ~exist('val','var') || isempty(val)
    val = 0;
  end
  if ~exist('wantcaution','var') || isempty(wantcaution)
    wantcaution = 1;
  end
else
  if isempty(val)
    val = 0;
  end
  if isempty(wantcaution)
    wantcaution = 1;
  end
end

% handle special case of y being scalar
if isscalar(y)
  if y==0
    f = repmat(val,size(x));
  else
    if wantcaution && abs(y) < 1e-5   % see allzero.m
      warning('abs value of divisor is less than 1e-5. we are treating the divisor as 0.');
      f = repmat(val,size(x));
    else
%REMOVED:
%      if abs(y) < 1e-5 && ~wantcaution
%        warning('abs value of divisor is less than 1e-5. we are treating the divisor as-is.');
%      end
      f = x./y;
    end
  end
else
  % do it
  bad = y==0;
  bad2 = abs(y) < 1e-5;  % see allzero.m
  if wantcaution && any(bad2(:) & ~bad(:))
    warning('abs value of one or more divisors is less than 1e-5. we are treating these divisors as 0.');
  end
%REMOVED:
%  if any(bad2 & ~bad) && ~wantcaution
%    warning('abs value of one or more divisors is less than 1e-5. we are treating the divisors as-is.');
%  end
  if wantcaution
    y(bad2) = 1;
    f = x./y;
    f(bad2) = val;
  else
    y(bad) = 1;
    f = x./y;
    f(bad) = val;
  end
end


function f = normalizemax(m,dim)

% function f = normalizemax(m,dim)
%
% <m> is a matrix
% <dim> (optional) is the dimension of <m> to operate upon.
%   default to 2 if <m> is a row vector and to 1 otherwise.
%   special case is 0 which means operate globally.
%
% divide <m> by the max value along some dimension (or globally).
%
% example:
% isequal(normalizemax([1 2 3]),[1/3 2/3 1])

% input
if ~exist('dim','var') || isempty(dim)
  dim = choose(isrowvector(m),2,1);
end

% do it
if dim==0
  f = m / max(m(:));
else
  f = bsxfun(@rdivide,m,max(m,[],dim));
end



function f = isrowvector(m)

% function f = isrowvector(m)
%
% <m> is a matrix
%
% return whether <m> is 1 x n where n >= 0.
% specifically:
%   f = isvector(m) & size(m,1)==1;
%
% example:
% isrowvector([1 2])
% isrowvector([1])
% isrowvector(zeros(1,0))
% ~isrowvector([])

f = isvector(m) & size(m,1)==1;



function f = choose(flag,yes,no)

% function f = choose(flag,yes,no)
%
% <flag> is a truth value (0 or 1)
% <yes> is something
% <no> is something
%
% if <flag>, return <yes>.  otherwise, return <no>.
%
% example:
% isequal(cellfun(@(x) choose(isempty(x),2,x),{[] 1}),[2 1])

if flag
  f = yes;
else
  f = no;
end



function f = squish(m,num)

% function f = squish(m,num)
%
% <m> is a matrix
% <num> is the positive number of initial dimensions to squish together
%
% return <m> squished.
%
% example:
% isequal(squish([1 2; 3 4],2),[1 3 2 4]')

% get the size of m
msize = [size(m) ones(1,num-ndims(m))];  % add ones to end if necessary

% calculate the new dimensions
newdim = [prod(msize(1:num)) msize(num+1:end)];

% do the reshape
f = reshape(m,[newdim 1]);  % tack on a 1 to handle the special case of squishing everything together



function f = mtimescell(m1,m2)

% function f = mtimescell(m1,m2)
%
% <m1> is A x B
% <m2> is a cell vector of matrices such that cat(1,m2{:}) is B x C
%
% simply return <m1>*cat(1,m2{:}) but do so in a way that doesn't cause 
% too much memory usage.
%
% example:
% x = randn(10,20);
% y = randn(20,200);
% result = x*y;
% result2 = mtimescell(x,splitmatrix(y,1,repmat(2,[1 10])));
% allzero(result-result2)

f = 0;
cnt = 0;
for q=1:length(m2)
  f = f + m1(:,cnt + (1:size(m2{q},1))) * m2{q};
  cnt = cnt + size(m2{q},1);
end



function m = catcell(dim,m)

% function m = catcell(dim,m)
%
% <dim> is the dimension to concatenate along
% <m> is a cell matrix
%
% simply return cat(dim,m{:}).  this function is useful because 
% MATLAB doesn't provide an easy way to apply "{:}" to an 
% arbitrary matrix.
%
% example:
% isequal(catcell(2,{1 2 3}),[1 2 3])

m = cat(dim,m{:});

% f = [];
% for p=1:numel(m)
%   if p == 1
%     f = m{p};
%   else
%     f = cat(dim,f,m{p});
%   end
%   m{p} = [];
% end
