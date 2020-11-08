function [N,crossValDelay,binCent] = runCrossValDelay(dOrig,p,verbose)

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
nBin = p.nBin;
histWidth = p.histWidth;
edges = [-histWidth/2:histWidth/nBin:histWidth/2];
binCent = edges(1:end-1)+2*pi/nBin/2;

for kInd = 1:length(unique(d.crossVal))
    [dTr{1,1},dTe{1,1},trInd{1,1},teInd{1,1},info] = splitTrainTest(d,crossVal(kInd));
    
    trDelayDiff(kInd,:) = wrapToPi( circ_mean(dTr{1,1}{2}.x(dTr{1,1}{2}.y==2,:),[],1) - circ_mean(dTr{1,1}{2}.x(dTr{1,1}{2}.y==1,:),[],1) );
    teDelayDiff(kInd,:) = wrapToPi( circ_mean(dTe{1,1}{2}.x(dTe{1,1}{2}.y==2,:),[],1) - circ_mean(dTe{1,1}{2}.x(dTe{1,1}{2}.y==1,:),[],1) );
    
    [N(kInd,:),~,bin] = histcounts(trDelayDiff(kInd,:),edges);
    binInd = nan(length(N(kInd,:)),length(bin));
    for i = 1:length(N(kInd,:))
        binInd(i,:) = bin==i;
        if isempty(find(binInd(i,:),1))
            crossValDelay(kInd,i) = nan;
        else
            crossValDelay(kInd,i) = circ_mean(teDelayDiff(kInd,logical(binInd(i,:))),[],2);
        end
        delay(kInd,i) = circ_mean(trDelayDiff(kInd,logical(binInd(i,:))),[],2);
    end
end

for i = 1:size(crossValDelay,2)
    tmp = crossValDelay(~isnan(crossValDelay(:,i)),i);
    if isempty(tmp)
        crossValDelay2(1,i) = nan;
    else
        crossValDelay2(1,i) = circ_mean(tmp);
    end
end
crossValDelay = crossValDelay2;
% crossValDelay = circ_mean(crossValDelay,[],1);
delay = circ_mean(delay,[],1);
N = mean(N,1);

if verbose
    figure('WindowStyle','docked');
    yyaxis left
    bar(binCent,N,1,'faceColor',get(gca,'YColor'),'edgeColor','w'); hold on
    xlim([edges(1) edges(end)])
    yyaxis right
    plot(binCent,crossValDelay,'-o'); hold on
    plot(binCent,delay,'-ro'); hold on
    plot(get(gca,'xlim'),[0 0],':')
end






function [dTr2,dTe2,trInd,teInd,info] = splitTrainTest(d,crossVal)

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
trPhase = angle(dTr.x);
tePhase = angle(dTe.x);
trMag = abs(dTr.x);
teMag = abs(dTe.x);
phaseShift = circ_mean([circ_mean(trPhase(dTr.y==1,:),[],1);circ_mean(trPhase(dTr.y==2,:),[],1)],[],1);
phaseShift_roi = repmat(circ_median(phaseShift,2),1,size(phaseShift,2));
magScale = mean(trMag,1);

% phaseShift = zeros(size(phaseShift));
% phaseShift_roi = zeros(size(phaseShift_roi));
% magScale = ones(size(magScale));

tr.magPol = dTr;
tr.magPol.x = trMag./repmat(magScale,size(trMag,1),1);
tr.delayPol = dTr;
tr.delayPol.x = wrapToPi(trPhase-repmat(phaseShift,size(trPhase,1),1));
tr.delayCartr = dTr; tr.delayCarti = dTr;
[tr.delayCartr.x,tr.delayCarti.x] = pol2cart(wrapToPi(trPhase-repmat(phaseShift,size(trPhase,1),1)),ones(size(trMag)));
tr.delayCartr_roiDelay = dTr; tr.delayCarti_roiDelay = dTr;
[tr.delayCartr_roiDelay.x,tr.delayCarti_roiDelay.x] = pol2cart(wrapToPi(trPhase-repmat(phaseShift_roi,size(trPhase,1),1)),ones(size(trMag)));
tr.cartr = dTr; tr.carti = dTr;
[tr.cartr.x,tr.carti.x] = pol2cart(wrapToPi(trPhase-repmat(phaseShift,size(trPhase,1),1)),trMag./repmat(magScale,size(trMag,1),1));
tr.cartr_roiDelay = dTr; tr.carti_roiDelay = dTr;
[tr.cartr_roiDelay.x,tr.carti_roiDelay.x] = pol2cart(wrapToPi(trPhase-repmat(phaseShift_roi,size(trPhase,1),1)),trMag./repmat(magScale,size(trMag,1),1));

te.magPol = dTe;
te.magPol.x = teMag./repmat(magScale,size(teMag,1),1);
te.delayPol = dTe;
te.delayPol.x = wrapToPi(tePhase-repmat(phaseShift,size(tePhase,1),1));
te.delayCartr = dTe; te.delayCarti = dTe;
[te.delayCartr.x,te.delayCarti.x] = pol2cart(wrapToPi(tePhase-repmat(phaseShift,size(tePhase,1),1)),ones(size(teMag)));
te.delayCartr_roiDelay = dTe; te.delayCarti_roiDelay = dTe;
[te.delayCartr_roiDelay.x,te.delayCarti_roiDelay.x] = pol2cart(wrapToPi(tePhase-repmat(phaseShift_roi,size(tePhase,1),1)),ones(size(teMag)));
te.cartr = dTe; te.carti = dTe;
[te.cartr.x,te.carti.x] = pol2cart(wrapToPi(tePhase-repmat(phaseShift,size(tePhase,1),1)),teMag./repmat(magScale,size(teMag,1),1));
te.cartr_roiDelay = dTe; te.carti_roiDelay = dTe;
[te.cartr_roiDelay.x,te.carti_roiDelay.x] = pol2cart(wrapToPi(tePhase-repmat(phaseShift_roi,size(tePhase,1),1)),teMag./repmat(magScale,size(teMag,1),1));

info = fields(tr)';
for infoInd = 1:length(info)
    dTr2{infoInd} = tr.(info{infoInd});
    dTr2{infoInd}.info = info{infoInd};
    dTr2{infoInd}.infoDim = 'vox';
    
    dTe2{infoInd} = te.(info{infoInd});
    dTe2{infoInd}.info = info{infoInd};
    dTe2{infoInd}.infoDim = 'vox';
end
clear tr te val


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
    %Z-score
    xTr = dTr{i}.x;
%     xShift = mean(xTr,1);
%     xTr = xTr-repmat(xShift,[size(xTr,1) 1]);
%     xScale = std(xTr,1);
%     xTr = xTr./repmat(xScale,[size(xTr,1) 1]);
    
    xTe = dTe{i}.x;
%     xTe = xTe-repmat(xShift,[size(xTe,1) 1]);
%     xTe = xTe./repmat(xScale,[size(xTe,1) 1]);
    
    %Train
    if length(p.C)>1
        dTr{i}.svmStruct = svmtrain(dTr{i}.y,xTr,['-s 0 -t 0 -c ' num2str(p.C(i)) ' -q']);
    else
        dTr{i}.svmStruct = svmtrain(dTr{i}.y,xTr,['-s 0 -t 0 -c ' num2str(p.C) ' -q']);
    end
%     dTr{i}.svmStruct = svmtrain(dTr{i}.y,xTr,'-s 0 -t 1 -d 2 -g 1 -c 1 -q');
%     dTr{i}.svmStruct = svmtrain(dTr{i}.y,xTr,'-s 0 -t 1 -d 2 -c 1 -q');
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
                        dIn{1}{infoInd}.acc_crossVal(subFoldInd,foldInd) = tmp(1);
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

