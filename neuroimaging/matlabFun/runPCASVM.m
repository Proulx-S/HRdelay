function [dTr4,dTe4] = runPCASVM(dOrig,p,verbose)

p.getPattern = true;
p.perm = p.doPerm;

d.x = dOrig.xData;
d.y = dOrig.label;
d.crossVal = dOrig.crossVal;
d.label = [1 2];
d.labelPairs = nchoosek(d.label,2);
d.getPattern = p.getPattern;
clear dOrig

%% Loop over k random folds
crossVal = unique(d.crossVal);
for kInd = 1:length(unique(d.crossVal))
    tic
    if verbose
        display(['Fold ' num2str(crossVal(kInd)) '/' num2str(length(crossVal))])
    end
    
    switch p.infoComb
        case {'pcaOnResp_lowLevel','pcaOnTS_lowLevel','pcaOnResp_pcaOnSpace'}
            curP = p;
            if isnumeric(p.C)
            else
                curP.C = 1;
            end
            
            %% Split train and test trials
            [dTr{1,1},dTe{1,1},trInd{1,1},teInd{1,1},info] = splitTrainTest(d,crossVal(kInd));
            
            %% SVM on low-level
            [dTr2{1,1},dTe2{1,1}] = trainNtest(dTr{1,1},dTe{1,1},p,0); clear dTr dTe
            
            %% Bring up one level
            [dTr3(kInd,1),dTe3(kInd,1)] = levelUp(d,dTr2,dTe2,trInd,teInd,info,p.infoComb); clear dTr2 dTe2
            
            %% SVM on higher level
            [dTr3{kInd,1},dTe3{kInd,1}] = trainNtest(dTr3{kInd,1},dTe3{kInd,1},p,0);
        case {'respSVM__pcaNone_pcaNone','respSVM__pcaRun_pcaNone','respSVM__pcaTime_pcaNone','respSVM__pcaNone_pcaSpace','respSVM__pcaRun_pcaSpace','respSVM__pcaTime_pcaSpace','respSVM__pcaNone_pcaRun','respSVM__pcaRun_pcaRun','respSVM__pcaTime_pcaRun'}
            curP = p;
            if isnumeric(p.C)
            else
                curP.C = 1;
            end
            
            %% Split train and test trials
            [dTr{1,1},dTe{1,1},trInd{1,1},teInd{1,1},info] = splitTrainTest(d,p,crossVal(kInd));
            
            %% SVM on low-level
            [dTr2{1,1},dTe2{1,1}] = trainNtest(dTr{1,1},dTe{1,1},p,0); clear dTr dTe
            
            %% Bring up one level
            [dTr3(kInd,1),dTe3(kInd,1)] = levelUp(d,dTr2,dTe2,trInd,teInd,info,p.infoComb); clear dTr2 dTe2
            
            %% SVM on higher level
            [dTr3{kInd,1},dTe3{kInd,1}] = trainNtest(dTr3{kInd,1},dTe3{kInd,1},p,0);
            
        case 'respPCA'
            %% Split train and test trials
            [dTr{1,1},dTe{1,1},trInd{1,1},teInd{1,1},info] = splitTrainTest(d,p,crossVal(kInd));
            
        otherwise
            error('XX')
    end
end
clear dTr dTe dTr2 dTe2

%% Rearange
for kInd = 1:size(dTr3,1)
    dTr4(kInd,1) = collapseSubFold(d,dTr3(kInd,:));
    dTe4(kInd,1) = collapseSubFold(d,dTe3(kInd,:));
    dTr3(kInd,:) = cell(size(dTr3{kInd,:}));
    dTe3(kInd,:) = cell(size(dTe3{kInd,:}));
end
clear dTr3 dTe3

dTr4 = collapseFold(dTr4);
dTe4 = collapseFold(dTe4);

dTr4 = replaceD(d,dTr4,0); dTr4 = dTr4{1}{1};
dTe4 = replaceD(d,dTe4); dTe4 = dTe4{1}{1};

if exist('cOpt','var')
    cOpt = collapseFoldcOpt(cOpt);
    dTr4.cOpt = cOpt{1};
end







function [dTr2,dTe2,trInd,teInd,info] = splitTrainTest(d,p,crossVal)

trInd = find(d.crossVal~=crossVal);
teInd = find(d.crossVal==crossVal);

dTr = rmfield(d,{'label'});
dTr.x = d.x(trInd,:,:);
dTr.y = d.y(trInd,:,:);
dTr.crossVal = d.crossVal(trInd);
dTe = rmfield(d,{'label'});
dTe.x = d.x(teInd,:,:);
dTe.y = d.y(teInd,:,:);
dTe.crossVal = d.crossVal(teInd);

%% Normalize phase and mag
shift = mean(dTr.x,1);
dTr.x = dTr.x - repmat(shift,[size(dTr.x,1) 1 1]);
dTe.x = dTe.x - repmat(shift,[size(dTe.x,1) 1 1]);
scale = std(dTr.x,[],1);
dTr.x = dTr.x ./ repmat(scale,[size(dTr.x,1) 1 1]);
dTe.x = dTe.x ./ repmat(scale,[size(dTe.x,1) 1 1]);

for infoInd = 1:size(dTr.x,3)
    if strcmp(p.infoComb,'respSVM__pcaNone_pcaNone')
        info{infoInd} = ['TR' num2str(infoInd)];
    else
        info{infoInd} = ['pc' num2str(infoInd)];
    end
    
    dTr2{infoInd} = dTr;
    dTr2{infoInd}.x = dTr.x(:,:,infoInd);
    dTr2{infoInd}.info = info{infoInd};
    dTr2{infoInd}.infoDim = 'vox';
    
    dTe2{infoInd} = dTe;
    dTe2{infoInd}.x = dTe.x(:,:,infoInd);
    dTe2{infoInd}.info = info{infoInd};
    dTe2{infoInd}.infoDim = 'vox';
end


function [dTr,dTe] = trainNtest(dTr,dTe,p,ALflag)
if ~exist('ALflag','var')
    ALflag = 0;
end
% %Reorganize data
% dTr_cur = dTr{1}; dTr_cur.x = [];
% dTe_cur = dTe{1}; dTe_cur.x = [];
% newxTr = [];
% newxTe = [];
for i = 1:length(dTr)
    xTr = dTr{i}.x;
    xTe = dTe{i}.x;
    
    %Train
    if length(p.C)>1
        dTr{i}.svmStruct = svmtrain(dTr{i}.y,xTr,['-s 0 -t 0 -c ' num2str(p.C(i)) ' -q']);
    else
        dTr{i}.svmStruct = svmtrain(dTr{i}.y,xTr,['-s 0 -t 0 -c ' num2str(p.C) ' -q']);
    end
    [dTr{i}.pred,acc,dTr{i}.yhat] = svmpredict(dTr{i}.y,xTr,dTr{i}.svmStruct,'-q');
    dTr{i}.acc = acc(1);
    
%     %%%%%%%%%%%%%%%%%%
%     [dTr{i}.w,dTr{i}.a] = getAandW(dTr{i}.svmStruct,dTr{i}.x);
%     
%      tmpxTr = mean(dTr{i}.x(dTr{i}.y==1,:),1)-mean(dTr{i}.x(dTr{i}.y==2,:),1);
%      tmpxTe = mean(dTe{i}.x(dTe{i}.y==1,:),1)-mean(dTe{i}.x(dTe{i}.y==2,:),1);
% %     tmpxTr = mean(dTr{i}.x(dTr{i}.y==2,:),1);
% %     tmpxTe = mean(dTe{i}.x(dTe{i}.y==2,:),1);
%     [dTr{i}.tmp.RHO_tr2te,dTr{i}.tmp.PVAL_tr2te] = corr(tmpxTr',tmpxTe');
%         [dTr{i}.tmp.RHO_tr2w,dTr{i}.tmp.PVAL_tr2w] = corr(tmpxTr',dTr{i}.w');
%         [dTr{i}.tmp.RHO_tr2a,dTr{i}.tmp.PVAL_tr2a] = corr(tmpxTr',dTr{i}.a');
%         [dTr{i}.tmp.RHO_te2w,dTr{i}.tmp.PVAL_te2w] = corr(tmpxTe',dTr{i}.w');
%         [dTr{i}.tmp.RHO_te2a,dTr{i}.tmp.PVAL_te2a] = corr(tmpxTe',dTr{i}.a');
%         [dTr{i}.tmp.RHO_te2w,dTr{i}.tmp.PVAL_te2w] = corr(tmpxTe',dTr{i}.w');
%         [dTr{i}.tmp.RHO_te2a,dTr{i}.tmp.PVAL_te2a] = corr(tmpxTe',dTr{i}.a');
% 
%         dTr{i}.tmp.N = length(tmpxTr);
%     %%%%%%%%%%%%%%%%%%
    
    
    %Test test data (for cross-validation)
    [dTe{i}.pred,acc,dTe{i}.yhat] = svmpredict(dTe{i}.y,xTe,dTr{i}.svmStruct,'-q');
    dTe{i}.acc = acc(1);
    
    %Revert anti-learning to learning
    if ALflag
        dTe{i}.yhat = dTe{i}.yhat*dTr{i}.ALrev;
        dTe{i}.pred = ((dTe{i}.pred-1.5)*dTr{i}.ALrev)+1.5;
        dTe{i}.acc = ((dTe{i}.acc-50)*dTr{i}.ALrev)+50;
    end
    
end

function [dTr2,dTe2] = catInfo(d,dTr,dTe,trInd,teInd,combInfo,infoStr)
if ~isempty(combInfo)
    for i  =1:length(dTr{1})
        allInfo{i} = dTr{1}{i}.info;
    end
    for i  =1:length(dTr)
        dTr2{i} = dTr{i}(ismember(allInfo,combInfo));
        dTe2{i} = dTe{i}(ismember(allInfo,combInfo));
    end
    dTr = dTr2; clear dTr2
    dTe = dTe2; clear dTe2
end

newTr = nan(size(d.y,1),length(dTr),length(dTr{1})*size(d.x,2));
newTe = nan(size(d.y,1),length(dTr),length(dTr{1})*size(d.x,2));
for kSubInd = 1:length(dTr)
    for infoInd = 1:length(dTr{1})
        curInd = (infoInd-1)*size(d.x,2)+1:(infoInd)*size(d.x,2);
        newTr(trInd{kSubInd},kSubInd,curInd) = permute(dTr{kSubInd}{infoInd}.x,[1 3 2]);
        newTe(teInd{kSubInd},kSubInd,curInd) = permute(dTe{kSubInd}{infoInd}.x,[1 3 2]);
    end
end

for kSubInd = 1:length(dTr)
    dTr2{kSubInd}{1} = dTr{kSubInd}{1};
    dTr2{kSubInd}{1}.svmStruct = [];
    dTr2{kSubInd}{1}.pred = [];
    dTr2{kSubInd}{1}.yhat = [];
    dTr2{kSubInd}{1}.acc = [];
    dTr2{kSubInd}{1}.ALrev = [];
    dTr2{kSubInd}{1}.x = squeeze(newTr(trInd{kSubInd},:,:)); % Not actually averaging across subFolds, rather compiling out the nans
    dTr2{kSubInd}{1}.info = infoStr;
    dTr2{kSubInd}{1}.infoDim = combInfo;
    dTr2{kSubInd}{1}.subLevel = dTr{kSubInd};
    

    dTe2{kSubInd}{1} = dTe{kSubInd}{1};
    dTe2{kSubInd}{1}.pred = [];
    dTe2{kSubInd}{1}.yhat = [];
    dTe2{kSubInd}{1}.acc = [];
    dTe2{kSubInd}{1}.x = squeeze(newTe(teInd{kSubInd},kSubInd,:));
    dTe2{kSubInd}{1}.info = infoStr;
    dTe2{kSubInd}{1}.infoDim = combInfo;
    dTe2{kSubInd}{1}.subLevel = dTe{kSubInd};
end

function [dTr2,dTe2] = levelUp(d,dTr,dTe,trInd,teInd,combInfo,infoStr,keepSubLevel)
if ~exist('keepSubLevel','var')
    keepSubLevel = 1;
end
if ~isempty(combInfo)
    for i  =1:length(dTr{1})
        allInfo{i} = dTr{1}{i}.info;
    end
    for i  =1:length(dTr)
        dTr2{i} = dTr{i}(ismember(allInfo,combInfo));
        dTe2{i} = dTe{i}(ismember(allInfo,combInfo));
    end
    dTr = dTr2; clear dTr2
    dTe = dTe2; clear dTe2
end

newTr = nan(size(d.y,1),length(dTr),length(dTr{1}));
newTe = nan(size(d.y,1),length(dTr),length(dTr{1}));

for kSubInd = 1:length(dTr)
    for infoInd = 1:length(dTr{1})
        newTr(trInd{kSubInd},kSubInd,infoInd) = dTr{kSubInd}{infoInd}.yhat;
        newTe(teInd{kSubInd},kSubInd,infoInd) = dTe{kSubInd}{infoInd}.yhat;
    end
end
for kSubInd = 1:length(dTr)
    dTr2{kSubInd}{1} = dTr{kSubInd}{1};
    dTr2{kSubInd}{1}.svmStruct = [];
    dTr2{kSubInd}{1}.pred = [];
    dTr2{kSubInd}{1}.yhat = [];
    dTr2{kSubInd}{1}.acc = [];
    dTr2{kSubInd}{1}.ALrev = [];
    dTr2{kSubInd}{1}.x = squeeze(nanmean(newTr(trInd{kSubInd},:,:),2)); % Not actually averaging across subFolds, rather compiling out the nans
    dTr2{kSubInd}{1}.info = infoStr;
    dTr2{kSubInd}{1}.infoDim = combInfo;
    if keepSubLevel
        dTr2{kSubInd}{1}.subLevel = dTr{kSubInd};
    elseif isfield(dTr2{kSubInd}{1},'subLevel')
        dTr2{kSubInd}{1} = rmfield(dTr2{kSubInd}{1},'subLevel');
    end
    
    dTe2{kSubInd}{1} = dTe{kSubInd}{1};
    dTe2{kSubInd}{1}.pred = [];
    dTe2{kSubInd}{1}.yhat = [];
    dTe2{kSubInd}{1}.acc = [];
    dTe2{kSubInd}{1}.x = squeeze(newTe(teInd{kSubInd},kSubInd,:));
    dTe2{kSubInd}{1}.info = infoStr;
    dTe2{kSubInd}{1}.infoDim = combInfo;
    if keepSubLevel
        dTe2{kSubInd}{1}.subLevel = dTe{kSubInd};
    elseif isfield(dTe2{kSubInd}{1},'subLevel')
        dTe2{kSubInd}{1} = rmfield(dTe2{kSubInd}{1},'subLevel');
    end
end

function [dTr2,dTe2] = voxelStplit(d,dTr,dTe,trInd,teInd,combInfo,infoStr)
for i  =1:length(dTr{1})
    allInfo{i} = dTr{1}{i}.info;
end
for i  =1:length(dTr)
        dTr2{i} = dTr{i}(ismember(allInfo,combInfo));
        dTe2{i} = dTe{i}(ismember(allInfo,combInfo));
end
dTr = dTr2; clear dTr2
dTe = dTe2; clear dTe2

for vox = 1:size(dTr{1}{1}.x,2)
    newTr = nan(size(d.y,1),length(dTr),length(dTr{1}));
    newTe = nan(size(d.y,1),length(dTr),length(dTr{1}));
    
    for kSubInd = 1:length(dTr)
        for infoInd = 1:length(dTr{1})
            newTr(trInd{kSubInd},kSubInd,infoInd) = dTr{kSubInd}{infoInd}.x(:,vox);
            newTe(teInd{kSubInd},kSubInd,infoInd) = dTe{kSubInd}{infoInd}.x(:,vox);
        end
    end
    
    for kSubInd = 1:length(dTr)
        dTr2{kSubInd}{vox} = dTr{kSubInd}{1};
        dTr2{kSubInd}{vox}.svmStruct = [];
        dTr2{kSubInd}{vox}.pred = [];
        dTr2{kSubInd}{vox}.yhat = [];
        dTr2{kSubInd}{vox}.acc = [];
        dTr2{kSubInd}{vox}.ALrev = [];
        dTr2{kSubInd}{vox}.x = squeeze(newTr(trInd{kSubInd},:,:)); % Not actually averaging across subFolds, rather compiling out the nans
        dTr2{kSubInd}{vox}.info = infoStr;
        dTr2{kSubInd}{vox}.infoDim = combInfo;
        dTr2{kSubInd}{vox}.subLevel = dTr{kSubInd};
        
        
        dTe2{kSubInd}{vox} = dTe{kSubInd}{1};
        dTe2{kSubInd}{vox}.pred = [];
        dTe2{kSubInd}{vox}.yhat = [];
        dTe2{kSubInd}{vox}.acc = [];
        dTe2{kSubInd}{vox}.x = squeeze(newTe(teInd{kSubInd},kSubInd,:));
        dTe2{kSubInd}{vox}.info = infoStr;
        dTe2{kSubInd}{vox}.infoDim = combInfo;
        dTe2{kSubInd}{vox}.subLevel = dTe{kSubInd};
    end
end



function dOut = collapseSubFold(d,dIn)
for kInd = 1:length(dIn)
    for infoInd = 1:length(dIn{kInd})
        allField = fields(dIn{kInd}{infoInd});
        for fieldInd = 1:length(allField)
            switch allField{fieldInd}
                case 'x'
%                     tmp = nan(size(d.(allField{fieldInd}),1),size(dIn{kInd}{infoInd}.(allField{fieldInd}),2));
%                     tmp(ismember(d.crossVal,dIn{kInd}{infoInd}.crossVal),:) = dIn{kInd}{infoInd}.(allField{fieldInd});
%                     dOut{infoInd}.(allField{fieldInd})(:,:,kInd) = tmp;
                case {'y','pred','yhat'}
                    tmp = nan(size(d.crossVal,1),size(dIn{kInd}{infoInd}.(allField{fieldInd}),2));
                    tmp(ismember(d.crossVal,dIn{kInd}{infoInd}.crossVal),:) = dIn{kInd}{infoInd}.(allField{fieldInd});
                    dOut{infoInd}.(allField{fieldInd})(:,kInd) = tmp;
                case {'d','crossVal','labelPairs','getPattern',}
                    dOut{infoInd}.(allField{fieldInd}) = d.(allField{fieldInd});
                case {'info','infoDim'}
                    dOut{infoInd}.(allField{fieldInd}) = dIn{kInd}{infoInd}.(allField{fieldInd});
                case {'svmStruct','acc'}
                    dOut{infoInd}.(allField{fieldInd})(kInd,1) = dIn{kInd}{infoInd}.(allField{fieldInd});
                case 'subLevel'
                    dOut{infoInd}.(allField{fieldInd}) = [];
                case {'w','a'}
                    dOut{infoInd}.(allField{fieldInd})(:,kInd) = dIn{kInd}{infoInd}.(allField{fieldInd})';
                case 'tmp'
                    allField2 = fields(dIn{kInd}{infoInd}.(allField{fieldInd}));
                    for i = 1:length(allField2)
                        dOut{infoInd}.(allField{fieldInd}).(allField2{i})(kInd,1) = dIn{kInd}{infoInd}.(allField{fieldInd}).(allField2{i});
                    end
                case 'ALrev'
                otherwise
                    error('XX')
            end
        end
    end
end
if isfield(dOut{infoInd},'subLevel')
    for infoInd = 1:length(dIn{kInd})
        for kInd = 1:length(dIn)
            dIn2{kInd} = dIn{kInd}{infoInd}.subLevel;
        end
        tmp = collapseSubFold(d,dIn2);
        dOut{infoInd}.subLevel = tmp{1};
    end
end

dOut = {dOut};


function dOut = collapseFold(dIn)

for kInd = 1:length(dIn)
    for infoInd = 1:length(dIn{kInd})
        allField = fields(dIn{kInd}{infoInd});
        for fieldInd = 1:length(allField)
            switch allField{fieldInd}
                case 'x'
                    dOut{infoInd}.(allField{fieldInd})(:,:,:,kInd) = dIn{kInd}{infoInd}.(allField{fieldInd});
                case {'y','pred','yhat'}
                    dOut{infoInd}.(allField{fieldInd})(:,:,kInd) = dIn{kInd}{infoInd}.(allField{fieldInd});
                case {'crossVal','labelPairs','getPattern',}
                    dOut{infoInd}.(allField{fieldInd}) = dIn{kInd}{infoInd}.(allField{fieldInd});
                case {'info','infoDim'}
                    dOut{infoInd}.(allField{fieldInd}) = dIn{kInd}{infoInd}.(allField{fieldInd});
                case {'svmStruct','acc','ALrev'}
                    dOut{infoInd}.(allField{fieldInd})(:,kInd) = dIn{kInd}{infoInd}.(allField{fieldInd});
                case 'subLevel'
                    dOut{infoInd}.(allField{fieldInd}) = [];
                case {'w','a'}
                    dOut{infoInd}.(allField{fieldInd})(:,:,kInd) = dIn{kInd}{infoInd}.(allField{fieldInd});
                case 'tmp'
                    allField2 = fields(dIn{kInd}{infoInd}.(allField{fieldInd}));
                    for i = 1:length(allField2)
                        dOut{infoInd}.(allField{fieldInd}).(allField2{i})(:,kInd) = dIn{kInd}{infoInd}.(allField{fieldInd}).(allField2{i});
                    end
                otherwise
                    error('XX')
            end
        end
    end
end
if isfield(dOut{infoInd},'subLevel')
    for infoInd = 1:length(dIn{kInd})
        for kInd = 1:length(dIn)
            dIn2{kInd} = dIn{kInd}{infoInd}.subLevel;
        end
        tmp = collapseFold(dIn2);
        dOut{infoInd}.subLevel = tmp{1};
    end
end
dOut = {dOut};





function dIn = replaceD(d,dIn,replace)
if ~exist('replace','var')
    replace = 1;
end

for infoInd = 1:length(dIn{1})
    allField = fields(dIn{1}{infoInd});
    for fieldInd = 1:length(allField)
        switch allField{fieldInd}
            case 'y'
                dIn{1}{infoInd}.(allField{fieldInd}) = nanmean(nanmean(dIn{1}{infoInd}.(allField{fieldInd}),3),2);
            case {'x','label','labelPairs','getPattern'}
                if replace
                    dIn{1}{infoInd}.(allField{fieldInd}) = d.(allField{fieldInd});
                else
                    dIn{1}{infoInd}.(allField{fieldInd}) = [];
                end
            case 'crossVal'
                dIn{1}{infoInd}.(allField{fieldInd}) = d.(allField{fieldInd});
            case {'info','infoDim','yhat','acc','ALrev'}
            case {'svmStruct','pred'} 
                dIn{1}{infoInd}.(allField{fieldInd}) = [];
            case 'subLevel'
                 tmp = replaceD(d,{dIn{1}{infoInd}.(allField{fieldInd})},0);
                 dIn{1}{infoInd}.(allField{fieldInd}) = tmp{1};
            otherwise
                error('XX')
        end
    end
    
    if ~isfield(dIn{1}{infoInd},'svmStruct')
        if isfield(dIn{1}{infoInd},'acc')
            for foldInd = 1:size(dIn{1}{infoInd}.acc,2)
                for subFoldInd = 1:size(dIn{1}{infoInd}.acc,1)
                    tmp = d.crossVal(~isnan(dIn{1}{infoInd}.yhat(:,subFoldInd,foldInd)));
                    if all(diff(tmp)==0)
                        if ~isempty(tmp)
                            dIn{1}{infoInd}.acc_crossVal(subFoldInd,foldInd) = tmp(1);
                        else
                            dIn{1}{infoInd}.acc_crossVal(subFoldInd,foldInd) = nan;
                        end
                    else
                        error('XX')
                    end
                end
            end
        end
    end
end

function [w,a] = getAandW(svmStruct,xTr)
%% Get w
% w
w = zeros(1,size(svmStruct.SVs,2));
for i = 1:size(svmStruct.SVs,1)
    w = w + svmStruct.sv_coef(i)*svmStruct.SVs(i,:);
end

%% Get a
% a
s = w*xTr';
% a = w*cov(xTr)/cov(s);
a = nan(size(w));
for dim = 1:size(xTr,2)
    tmp = cov(xTr(:,dim),s);
    a(dim) = tmp(1,2);
end


function [dTr3,dTe3] = doInfoComb(d,p,crossVal,cOpt,lowLevelCflag)
if ~exist('lowLevelCflag','var')
    lowLevelCflag = 0;
end
%% Split train and test trials
[dTr{1,1},dTe{1,1},trInd{1,1},teInd{1,1},info] = splitTrainTest(d,crossVal);

%% SVM on low-level
curP = p;
%adjust C if needed
if exist('cOpt','var') && ~isempty(cOpt)
    if lowLevelCflag
        curP.C = [];
        for i = 1:length(cOpt)
            curP.C = [curP.C cOpt{i}.subLevel{1}.c1 cOpt{i}.subLevel{2}.c1];
        end
    else
        curP.C = [];
        for i = 1:length(cOpt)
            curP.C = [curP.C repmat(cOpt{i}.c1,[1 2])];
        end
    end
else
    if isnumeric(p.C)
    else
        curP.C = p.C.Cmax;
    end
end
[dTr{1,1},dTe{1,1}] = trainNtest(dTr{1,1},dTe{1,1},curP,0);

%% Collapse info
info = reshape(info,[2 length(info)/2])';
infoComb = {'pol' 'delayCart' 'delayCart_roiDelay' 'cart' 'cart_roiDelay'}';
for infoInd = 1:size(info,1)
    [dTr_tmp,dTe_tmp] = levelUp(d,dTr(1,:),dTe(1,:),trInd(1,:),teInd(1,:),info(infoInd,:),infoComb{infoInd});
    dTr2{1,1}(infoInd) = dTr_tmp{1,1};
    dTe2{1,1}(infoInd) = dTe_tmp{1,1};
end
clear dTr dTe dTr_tmp dTe_tmp

if lowLevelCflag
    dTr3 = dTr2;
    dTe3 = dTe2;
else
    %adjust C if needed
    if exist('cOpt','var') && ~isempty(cOpt)
        curP.C = [];
        for i = 1:length(cOpt)
            curP.C(i) = cOpt{i}.c2;
        end
    else
        curP.C = 1;
    end
    
    [dTr2{1,1},dTe2{1,1}] = trainNtest(dTr2{1,1},dTe2{1,1},curP,0);
    
    %% Bring up one level
    [dTr3(1,1),dTe3(1,1)] = levelUp(d,dTr2,dTe2,trInd,teInd,[],p.infoComb); clear dTr2 dTe2
    
    %% SVM on higher level
    curP.C = 1;
    [dTr3{1,1},dTe3{1,1}] = trainNtest(dTr3{1,1},dTe3{1,1},curP,0);
end

function [dTr3,dTe3] = infoComb_Copt(d,p,crossVal,c1,c2,lastLevel)
if ~exist('lastLevel','var')
    lastLevel = 0;
end
%% Split train and test trials
[dTr{1,1},dTe{1,1},trInd{1,1},teInd{1,1},info] = splitTrainTest(d,crossVal);

%% SVM on low-level
curP = p;
curP.C = c1;
[dTr{1,1},dTe{1,1}] = trainNtest(dTr{1,1},dTe{1,1},curP,0);

%% Collapse info
info = reshape(info,[2 length(info)/2])';
infoComb = {'pol' 'delayCart' 'delayCart_roiDelay' 'cart' 'cart_roiDelay'}';
for infoInd = 1:size(info,1)
    [dTr_tmp,dTe_tmp] = levelUp(d,dTr(1,:),dTe(1,:),trInd(1,:),teInd(1,:),info(infoInd,:),infoComb{infoInd});
    dTr2{1,1}(infoInd) = dTr_tmp{1,1};
    dTe2{1,1}(infoInd) = dTe_tmp{1,1};
end
clear dTr dTe dTr_tmp dTe_tmp
for cInd = 1:length(c2)
    curP.C = c2(cInd);
    [dTr3{1,cInd},dTe3{1,cInd}] = trainNtest(dTr2{1,1},dTe2{1,1},curP,0);
end

% %% Bring up one level
% [dTr3,dTe3] = levelUp(d,dTr2,dTe2,trInd,teInd,[],p.infoComb); clear dTr2 dTe2
% 
% %% SVM on higher level
% % p.cList2 = p.cList;
% for cInd = 1:length(p.cList2)
%     curP.C = p.cList2(cInd);
%     [dTr4{1,cInd},dTe4{1,cInd}] = trainNtest(dTr3{1,1},dTe3{1,1},curP,0);
% end
% clear dTr3 dTe3

function cOpt = infoComb_Copt__sub(d,p,crossVal,kInd,c1,c2)
toRemove = d.crossVal==crossVal(kInd);
d.x(toRemove,:) = [];
d.y(toRemove,:) = [];
d.crossVal(toRemove,:) = [];
crossVal(kInd) = [];
tmpCrossVal = nan(size(d.crossVal));
for i = 1:length(crossVal)
    tmpCrossVal(d.crossVal==crossVal(i))=i;
end
d.crossVal = tmpCrossVal; clear tmpCrossVal

for kInd = 1:length(crossVal)
    [dTr3(kInd,:),dTe3(kInd,:)] = infoComb_Copt(d,p,kInd,c1,c2);
end
%% Rearange
for cInd = 1:size(dTr3,2)
    for kInd = 1:size(dTr3,1)
        dTr4(kInd,cInd) = collapseSubFold(d,dTr3(kInd,cInd));
        dTe4(kInd,cInd) = collapseSubFold(d,dTe3(kInd,cInd));
        dTr3(kInd,cInd) = cell(1,1);
        dTe3(kInd,cInd) = cell(1,1);
    end
end
clear dTr3 dTe3

for cInd = 1:size(dTr4,2)
    dTr5(1,cInd) = collapseFold(dTr4(:,cInd));
    dTe5(1,cInd) = collapseFold(dTe4(:,cInd));
    dTr4(:,cInd) = cell(size(dTr4(:,cInd)));
    dTe4(:,cInd) = cell(size(dTe4(:,cInd)));
end
clear dTr4 dTe4

for cInd = 1:size(dTr5,2)
    cOpt{1,cInd} = getDistT(dTe5{1,cInd});
end

function out = getDistT(dTe)
for anaInd = 1:length(dTe)
    yhat = nanmean(dTe{anaInd}.yhat,3);
    [~,~,~,STATS] = ttest2(yhat(1:end/2),yhat(end/2+1:end));
    out{anaInd}.distT = STATS.tstat;
    out{anaInd}.info = dTe{anaInd}.info;
    if isfield(dTe{anaInd},'subLevel')
        out{anaInd}.subLevel = getDistT(dTe{anaInd}.subLevel);
    end
end


function cOpt2 = compileCopt(cOpt,c1,c2)
%Compile
for anaInd = 1:size(cOpt{1,1},2)
    for i = 1:size(cOpt,1)
        for ii = 1:size(cOpt,2)
            cOpt2{anaInd}.distT(i,ii) = cOpt{i,ii}{anaInd}.distT;
            cOpt2{anaInd}.info = cOpt{i,ii}{anaInd}.info;
            for subLevelInd = 1:length(cOpt{i,ii}{anaInd}.subLevel)
                cOpt2{anaInd}.subLevel{subLevelInd}.distT(i,ii) = cOpt{i,ii}{anaInd}.subLevel{subLevelInd}.distT;
                cOpt2{anaInd}.subLevel{subLevelInd}.info = cOpt{i,ii}{anaInd}.subLevel{subLevelInd}.info;
                if isfield(cOpt{i,ii}{anaInd}.subLevel{subLevelInd},'subLevel')
                    for subLevelInd2 = 1:length(cOpt{i,ii}{anaInd}.subLevel{subLevelInd}.subLevel)
                        cOpt2{anaInd}.subLevel{subLevelInd}.subLevel{subLevelInd2}.distT(i,ii) = cOpt{i,ii}{anaInd}.subLevel{subLevelInd}.subLevel{subLevelInd2}.distT;
                        cOpt2{anaInd}.subLevel{subLevelInd}.subLevel{subLevelInd2}.info = cOpt{i,ii}{anaInd}.subLevel{subLevelInd}.subLevel{subLevelInd2}.info;
                    end
                end
            end
        end
    end
    
    %Extract max
    tmp = cOpt2{anaInd}.distT;
    tmp = tmp==max(max(tmp)); [I,J] = find(tmp); [~,b] = min(abs(1-sqrt(c1(I).^2+c2(J).^2)));
    cOpt2{anaInd}.c1 = c1(I(b));
    cOpt2{anaInd}.c2 = c2(J(b));
    for subLevelInd = 1:length(cOpt{i,ii}{anaInd}.subLevel)
        tmp = cOpt2{anaInd}.subLevel{subLevelInd}.distT(:,1);
        tmp = tmp==max(tmp);
        [~,b] = min(abs(1-c1(tmp)));
        cOpt2{anaInd}.subLevel{subLevelInd}.c1 = c1(b);
    end
end


