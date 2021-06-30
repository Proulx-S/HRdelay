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
    mtimescell(olsmatrix2(design), ...
    cellfun(@(x) squish(x,3)',X,'UniformOutput',0));  % regressors x voxels
res.betas = permute(reshape(res.betas,[size(res.betas,1) xyzSz]),[2 3 4 1]);
res.design = designUncensored;
res.censorPts = censorPts;
res.designInfo = designInfo;







% % OLS Reduced-Model 1
% modelInd = ~ismember(designInfo(2,:),{'cond1' 'cond2' 'cond3'});
% dataInd = ~catcell(1,d.censorPts);
% display('computing OLS')
% betas = ...
%     mtimescell(olsmatrix2(design(dataInd,modelInd)), ...
%     cellfun(@(x) squish(x(:,:,:,~d.censorPts{1}),3)',d.data,'UniformOutput',0));  % regressors x voxels
% 
% betas = permute(reshape(betas,[size(design(:,modelInd),2) p.xyzSz]),[2 3 4 1]);
% 
% ind = ismember(designInfo(4,modelInd),{'poly0'});
% base = betas(:,:,:,ind);



% 
% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MIXED-EFFECT FIT
% tic
% fprintf('*** MIXED(random)-EFFECT FIT ***\n');
% %% Construct design matrix
% % Polynomials and motion
% poly = cell(1,numruns);
% polyNmotion = cell(1,numruns);
% constant = cell(1,numruns);
% 
% numtime = size(data{p},dimtime)*splitIn;
% for p=1:numruns/splitIn
%     tmp = constructpolynomialmatrix(numtime,0:opt.maxpolydeg(p),0);
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
% %Reformat to remerge motion regressors that were split before
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
% % Conditions
% condLabel = cat(1,ones(numruns/3,1)*1,ones(numruns/3,1)*2,ones(numruns/3,1)*3);
% sessLabel = cell2mat(opt.sessionLabel)';
% 
% condDesign1 = cell(1,3);
% polyNmotionDesign1 = cell(1,3);
% constantDesign1 = cell(1,3);
% condDesign2 = cell(1,3);
% polyNmotionDesign2 = cell(1,3);
% constantDesign2 = cell(1,3);
% for cond = 1:3
%     sess=1;
%     condInd = condLabel==cond & sessLabel==sess;
%     condInd_polyNmotion = condInd(1:opt.splitedIn:end);
%     condDesign1{cond} = blkdiag(cache.rawdesign{condInd});
%     
%     polyNmotionDesign1{cond} = blkdiag(polyNmotion{condInd_polyNmotion});
%     constantDesign1{cond} = blkdiag(constant{condInd_polyNmotion});
%     
%     
%     sess=2;
%     condInd = condLabel==cond & sessLabel==sess;
%     condInd_polyNmotion = condInd(1:opt.splitedIn:end);
%     condDesign2{cond} = blkdiag(cache.rawdesign{condInd});
%     
%     polyNmotionDesign2{cond} = blkdiag(polyNmotion{condInd_polyNmotion});
%     constantDesign2{cond} = blkdiag(constant{condInd_polyNmotion});
% end
% condDesign = reshape(cat(1,condDesign1,condDesign2),1,6);
% polyNmotionDesign = reshape(cat(1,polyNmotionDesign1,polyNmotionDesign2),1,6);
% constantDesign = reshape(cat(1,constantDesign1,constantDesign2),1,6);
% 
% 
% % Merge conditions, polynomials and motion
% condDesign = blkdiag(condDesign{:});
% % sessDesign = catcell(1,sessDesign); sessDesign(:,1:2) = [];
% polyNmotionDesign = blkdiag(polyNmotionDesign{:});
% constantDesign = blkdiag(constantDesign{:});
% 
% results.OLS.mixed.designmatrix = cat(2,condDesign,polyNmotionDesign);% time x regressors
% % results.OLS.mixed.designmatrixCond3 = cat(2,condDesign(end*2/3+1:end,:),polyNmotionDesign(end*2/3+1:end,:));% time x regressors
% 
% % Extract some usefull params
% results.OLS.mixed.designmatrixPieces.cond = cat(2,condDesign,zeros(size(polyNmotionDesign)));
% results.OLS.mixed.designmatrixPieces.motion = cat(2,zeros(size(condDesign)),polyNmotionDesign);
% results.OLS.mixed.designmatrixPieces.constant = cat(2,zeros(size(condDesign)),constantDesign);
% %deal with higher order poly
% ind = find(any(results.OLS.mixed.designmatrixPieces.constant,1));
% for i = 1:unique(opt.maxpolydeg)
%     results.OLS.mixed.designmatrixPieces.(['p' num2str(i)]) = zeros(size(results.OLS.mixed.designmatrixPieces.constant));
%     results.OLS.mixed.designmatrixPieces.(['p' num2str(i)])(:,ind+i) = results.OLS.mixed.designmatrix(:,ind+i);
% end
% 
% %     close all
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.mixed.designmatrix)
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.mixed.designmatrixPieces.cond)
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.mixed.designmatrixPieces.motion)
% %     figure('WindowStyle','docked'); colormap gray
% %     imagesc(results.OLS.mixed.designmatrixPieces.constant,[-1 1])
% 
% 
% 
% %remove constant
% ind = any(results.OLS.mixed.designmatrixPieces.constant,1);
% results.OLS.mixed.designmatrix(:,ind) = [];
% fieldList = fields(results.OLS.mixed.designmatrixPieces);
% for i = 1:length(fieldList)
%     results.OLS.mixed.designmatrixPieces.(fieldList{i})(:,ind) = [];
% end
% results.OLS.mixed.designmatrixPieces = rmfield(results.OLS.mixed.designmatrixPieces,'constant');
% 
% %% Estimate parameters with OLS
% display('computing OLS')
% results.OLS.mixed.parameters = ...
%     mtimescell(olsmatrix2(results.OLS.mixed.designmatrix), ...
%     cellfun(@(x) squish(x,dimdata)',data,'UniformOutput',0));  % regressors x voxels
% 
% for i = 1:length(opt.sessionLabel)
%     respTmp(:,:,i) = results.OLS.mixed.parameters((i-1)*12+1:i*12,:);
% end
% for i = 1:unique(opt.maxpolydeg)
%     ind = any(results.OLS.mixed.designmatrixPieces.(['p' num2str(i)]),1);
%     polyPar.(['p' num2str(i)]) = permute(results.OLS.mixed.parameters(ind,:),[3 2 1]);
% end
% results.OLS.mixed = rmfield(results.OLS.mixed,'parameters');
% 
% respMixed1 = respTmp(:,:,cell2mat(opt.sessionLabel)==1);
% respMixed1 = cat(4,respMixed1(:,:,1:end/3),respMixed1(:,:,end/3+1:end/3*2),respMixed1(:,:,end/3*2+1:end));
% for i = 1:unique(opt.maxpolydeg)
%     polyPar1.(['p' num2str(i)]) = polyPar.(['p' num2str(i)])(:,:,cell2mat(opt.sessionLabel)==1);
%     polyPar1.(['p' num2str(i)]) = cat(4,polyPar1.(['p' num2str(i)])(:,:,1:end/3),polyPar1.(['p' num2str(i)])(:,:,end/3+1:end/3*2),polyPar1.(['p' num2str(i)])(:,:,end/3*2+1:end));
% end
% respMixed2 = respTmp(:,:,cell2mat(opt.sessionLabel)==2);
% respMixed2 = cat(4,respMixed2(:,:,1:end/3),respMixed2(:,:,end/3+1:end/3*2),respMixed2(:,:,end/3*2+1:end));
% for i = 1:unique(opt.maxpolydeg)
%     polyPar2.(['p' num2str(i)]) = polyPar.(['p' num2str(i)])(:,:,cell2mat(opt.sessionLabel)==2);
%     polyPar2.(['p' num2str(i)]) = cat(4,polyPar2.(['p' num2str(i)])(:,:,1:end/3),polyPar2.(['p' num2str(i)])(:,:,end/3+1:end/3*2),polyPar2.(['p' num2str(i)])(:,:,end/3*2+1:end));
% end
% clear respTmp polyPar
% 
% % %Percent BOLD
% % respMixed1 = (respMixed1 - repmat(mean(respMixed1,1),[12 1 1 1]))./repmat(mean(respMixed1,1),[12 1 1 1]);
% % respMixed2 = (respMixed2 - repmat(mean(respMixed2,1),[12 1 1 1]))./repmat(mean(respMixed2,1),[12 1 1 1]);
% 
% for run = 1:size(respMixed1,3)
%     for cond = 1:size(respMixed1,4)
%         results.OLS.mixed.sess1.resp(:,:,:,:,run,cond) =   reshape(respMixed1(:,:,run,cond)',  [xyzsize size(respMixed1(:,:,run,cond),1)]);
%         for i = 1:unique(opt.maxpolydeg)
%             results.OLS.mixed.sess1.(['p' num2str(i)])(:,:,:,:,run,cond) =   reshape(polyPar1.(['p' num2str(i)])(:,:,run,cond)',  [xyzsize size(polyPar1.(['p' num2str(i)])(:,:,run,cond),1)]);
%         end
%     end
% end
% for run = 1:size(respMixed2,3)
%     for cond = 1:size(respMixed2,4)
%         results.OLS.mixed.sess2.resp(:,:,:,:,run,cond) =   reshape(respMixed2(:,:,run,cond)',  [xyzsize size(respMixed2(:,:,run,cond),1)]);
%         for i = 1:unique(opt.maxpolydeg)
%             results.OLS.mixed.sess2.(['p' num2str(i)])(:,:,:,:,run,cond) =   reshape(polyPar2.(['p' num2str(i)])(:,:,run,cond)',  [xyzsize size(polyPar2.(['p' num2str(i)])(:,:,run,cond),1)]);
%         end
%     end
% end
% 
