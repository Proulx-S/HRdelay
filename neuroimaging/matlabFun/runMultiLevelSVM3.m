function [rTr,rTe] = runMultiLevelSVM3(dOrig,p)
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
    display(['Fold ' num2str(crossVal(kInd)) '/' num2str(length(crossVal))])
    
    %% Split train and test trials
    teInd = find(d.crossVal==crossVal(kInd));
    trInd = find(d.crossVal~=crossVal(kInd));
    
    dTr = rmfield(d,{'label'});
    dTr.x = d.x(trInd,:,:);
    dTr.y = d.y(trInd,:,:);
    dTr.crossVal = d.crossVal(trInd);
    dTe = rmfield(d,{'label'});
    dTe.x = d.x(teInd,:,:);
    dTe.y = d.y(teInd,:,:);
    dTe.crossVal = d.crossVal(teInd);
    
    %% Normalize phase and mag
    trPhase = angle(dTr.x); tePhase = angle(dTe.x);
    trMag = abs(dTr.x); teMag = abs(dTe.x);
    phaseShift = circ_mean([circ_mean(trPhase(dTr.y==1,:),[],1);circ_mean(trPhase(dTr.y==2,:),[],1)],[],1);
    phaseShift_roi = repmat(circ_median(phaseShift,2),1,size(phaseShift,2));
    magScale = mean(trMag,1);
    
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
    
    %% Run multilevel
    switch p.infoComb
        case {'ALrev1'}
            %Store voxel vector for each information type
            info = fields(tr)';
            for infoInd = 1:length(info)
                dTr2{infoInd} = tr.(info{infoInd});
                dTr2{infoInd}.info = info{infoInd};
                dTr2{infoInd}.infoDim = 'vox';
                
                dTe2{infoInd} = te.(info{infoInd});
                dTe2{infoInd}.info = info{infoInd};
                dTe2{infoInd}.infoDim = 'vox';
            end
            clear tr te
            
            % Anti-learning revert
            curP = p;
            curP.subCrossVal = [];
            curP.subCrossVal.train = 1;
            curP.subCrossVal.test = 1;
            
            for infoInd = 1:length(info)
                curInfo = [info{infoInd} '__' curP.infoComb];
                combInfo = info(infoInd);
                combInfoInd = ismember(info,combInfo);
                [dTr2{end+1}, dTe2{end+1}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),curP,combInfo,curInfo,'vox');
                info{end+1} = curInfo;
            end
            
            
            
            
            
            % Collapse all info
            curP = p;
            curP.subCrossVal = [];
            curP.subCrossVal.train = 1;
            curP.subCrossVal.test = 1;
            
            curInfo = curP.infoComb;
            combInfo = info;
            combInfoInd = ismember(info,combInfo);
            [dTr3{kInd}, dTe3{kInd}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),curP,combInfo,curInfo,'infoOrVox');
            
            
            curP = p;
            curP.subCrossVal = [];
            curP.subCrossVal.train = 1;
            curP.subCrossVal.test = 1;
            
            dTr3{kInd} = trainSVMatLevel2(dTr3{kInd},curP);
            dTe3{kInd} = testSVMatLevel2(dTe3{kInd},dTr3{kInd});
            clear dTr2 dTe2
            
        otherwise
            error('X')
    end
    toc
    
end

%% Remove unnecessary data
for kInd = 1:length(crossVal)
    dTr3{kInd} = removex(dTr3{kInd});
    dTe3{kInd} = removex(dTe3{kInd});
end

%% Compile folds
keyboard
%Need to make sure the subFold yhat are compiled appropriately for te
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rTr = compileFolds(dTr3,d);
rTe = compileFolds(dTe3,d);


%% Optionally explore anti-learning
exploreAL3(rTr,'train')
exploreAL3(rTe,'test')
exploreAL4(rTr,rTe,[3 7])
exploreAL5(rTr,rTe,[3 7])


rTr = averageFolds(rTr);
rTe = averageFolds(rTe);




function r = compileFolds(dFold,d)
r = dFold{1};
r.y = d.y;
r.crossVal = d.crossVal;
if isfield(r,'bias')
    r = rmfield(r,'bias');
end
r.pred = nan(size(d.y,1),length(dFold));
r.yhat = nan(size(d.y,1),length(dFold));
if isfield(r,'acc')
    r.acc = nan(1,length(dFold));
end
if isfield(r,'w')
    r.w = nan(length(dFold),size(r.w,2));
    r.a = nan(length(dFold),size(r.a,2));
end

for k = 1:length(dFold)
    if isfield(r,'norm')
        r.norm(:,k) = dFold{k}.norm;
    end
    if isfield(r,'svmStruct')
        r.svmStruct(:,k) = dFold{k}.svmStruct;
    end
    ind = ismember(d.crossVal,dFold{k}.crossVal);
    r.pred(ind,k) = dFold{k}.pred;
    r.yhat(ind,k) = dFold{k}.yhat;
    if isfield(r,'acc')
        r.acc(1,k) = dFold{k}.acc;
    end
    if isfield(r,'w')
        r.w(k,:) = dFold{k}.w;
        r.a(k,:) = dFold{k}.a;
    end
end
% if isfield(r,'acc')
%     r.acc = mean(r.acc,2);
% %     r.acc = sum(r.acc)./length(r.y);
%     %     r.acc = r.acc./(length(r.y)/length(dFold));
% end
% if isfield(r,'w')
%     r.w= mean(r.w,1);
%     r.a= mean(r.a,1);
% end
% r.yhat = nanmean(r.yhat,2);
% r.pred = nanmean(r.pred,2);

if isfield(r,'subLevel')
    for subInd = 1:length(dFold{1}.subLevel)
        dFold2 = cell(1,length(dFold));
        for k = 1:length(dFold)
            dFold2{k} = dFold{k}.subLevel{subInd};
        end
        r.subLevel{subInd} = compileFolds(dFold2,d);
    end
end

function r = averageFolds(r)
if isfield(r,'acc')
    r.acc = mean(r.acc,2);
end
if isfield(r,'w')
    r.w= mean(r.w,1);
    r.a= mean(r.a,1);
end
r.yhat = nanmean(r.yhat,2);
r.pred = nanmean(r.pred,2);

if isfield(r,'subLevel')
    for subInd = 1:length(r.subLevel)
        r.subLevel{subInd} = averageFolds(r.subLevel{subInd});
    end
end


function d = removex(d)
if isfield(d,'subLevel')
    for i = 1:length(d.subLevel)
        d.subLevel{i} = removex(d.subLevel{i});
    end
end
d.x = [];

function [dTr_cur,dTe_cur] = svmReduction(dTr,dTe,p,combFields,curInfo,redDim)

switch redDim
    case 'info'
        error('double-check that')
        %Reorganize data
        dTr_cur = dTr{1}; dTr_cur.x = [];
        dTe_cur = dTe{1}; dTe_cur.x = [];
        xTr = [];
        xTe = [];
        for i = 1:length(dTr)
            xTr = cat(3,xTr,dTr{i}.x);
            xTe = cat(3,xTe,dTe{i}.x);
        end
        xTr = permute(xTr,[1 3 2]);
        xTe = permute(xTe,[1 3 2]);
        newxTr = [];
        newxTe = [];
        for i = 1:size(xTr,3)
            %Train
            dTr_cur_cur{i} = dTr_cur;
            dTr_cur_cur{i}.x = xTr(:,:,i);
            dTr_cur_cur{i} = trainSVMatLevel2(dTr_cur_cur{i},p);
            dTr_cur_cur{i}.info = '@vox';
            dTr_cur_cur{i}.infoDim = combFields;
            newxTr = cat(2,newxTr,dTr_cur_cur{i}.yhat);
            %Test
            dTe_cur_cur{i} = dTe_cur;
            dTe_cur_cur{i}.x = xTe(:,:,i);
            dTe_cur_cur{i} = testSVMatLevel2(dTe_cur_cur{i},dTr_cur_cur{i});
            dTe_cur_cur{i}.info = '@vox';
            dTe_cur_cur{i}.infoDim = combFields;
            newxTe = cat(2,newxTe,nanmean(dTe_cur_cur{i}.yhat,2));
        end
        %Store data
        dTr_cur.x = newxTr;
        dTr_cur.info = curInfo;
        dTr_cur.infoDim = 'vox';
        dTr_cur.subLevel = dTr_cur_cur;
        dTe_cur.x = newxTe;
        dTe_cur.info = curInfo;
        dTe_cur.infoDim = 'vox';
        dTe_cur.subLevel = dTe_cur_cur;
        
    case 'vox'
        %Reorganize data
        dTr_cur = dTr{1}; dTr_cur.x = [];
        dTe_cur = dTe{1}; dTe_cur.x = [];
        xTr = [];
        xTe = [];
        for i = 1:length(dTr)
            xTr = cat(3,xTr,dTr{i}.x);
            xTe = cat(3,xTe,dTe{i}.x);
        end
        newxTr = [];
        newxTe = [];
        for i = 1:size(xTr,3)
            %Train
            dTr_cur_cur{i} = dTr{i};
            dTr_cur_cur{i}.x = xTr(:,:,i);
            dTr_cur_cur{i} = trainSVMatLevel2(dTr_cur_cur{i},p);
            newxTr = cat(2,newxTr,dTr_cur_cur{i}.yhat);
            %Test
            dTe_cur_cur{i} = dTe{i};
            dTe_cur_cur{i}.x = xTe(:,:,i);
            dTe_cur_cur{i} = testSVMatLevel2(dTe_cur_cur{i},dTr_cur_cur{i});
            newxTe = cat(2,newxTe,nanmean(dTe_cur_cur{i}.yhat,2));
        end
        %Store data
        dTr_cur.x = newxTr;
        dTr_cur.info = curInfo;
        dTr_cur.infoDim = combFields;
        dTr_cur.subLevel = dTr_cur_cur;
        dTe_cur.x = newxTe;
        dTe_cur.info = curInfo;
        dTe_cur.infoDim = combFields;
        dTe_cur.subLevel = dTe_cur_cur;
        
    case {'infoOrVox','voxOrInfo'}
        %Reorganize data
        dTr_cur = dTr{1}; dTr_cur.x = [];
        dTe_cur = dTe{1}; dTe_cur.x = [];
        newxTr = [];
        newxTe = [];
        for i = 1:length(dTr)
            %Train
            dTr{i} = trainSVMatLevel2(dTr{i},p);
            newxTr = cat(2,newxTr,dTr{i}.yhat);
            %Test
            dTe{i} = testSVMatLevel2(dTe{i},dTr{i});
            newxTe = cat(2,newxTe,nanmean(dTe{i}.yhat,2));
        end
        %Store data
        dTr_cur.x = newxTr;
        dTr_cur.info = curInfo;
        dTr_cur.infoDim = combFields;
        dTr_cur.subLevel = dTr;
        dTe_cur.x = newxTe;
        dTe_cur.info = curInfo;
        dTe_cur.infoDim = combFields;
        dTe_cur.subLevel = dTe;
        
    otherwise
        error('X')
end
%
%
% function [dTr_cur,dTe_cur] = svmReduction2(dTr,dTe,p,combFields,curInfo,redDim)
%
% switch redDim
%     case 'info'
%         error('double-check that')
%         %Reorganize data
%         dTr_cur = dTr{1}; dTr_cur.x = [];
%         dTe_cur = dTe{1}; dTe_cur.x = [];
%         xTr = [];
%         xTe = [];
%         for i = 1:length(dTr)
%             xTr = cat(3,xTr,dTr{i}.x);
%             xTe = cat(3,xTe,dTe{i}.x);
%         end
%         xTr = permute(xTr,[1 3 2]);
%         xTe = permute(xTe,[1 3 2]);
%         newxTr = [];
%         newxTe = [];
%         for i = 1:size(xTr,3)
%             %Train
%             dTr_cur_cur{i} = dTr_cur;
%             dTr_cur_cur{i}.x = xTr(:,:,i);
%             dTr_cur_cur{i} = trainSVMatLevel(dTr_cur_cur{i},p);
%             dTr_cur_cur{i}.info = '@vox';
%             dTr_cur_cur{i}.infoDim = combFields;
%             newxTr = cat(2,newxTr,dTr_cur_cur{i}.yhat);
%             %Test
%             dTe_cur_cur{i} = dTe_cur;
%             dTe_cur_cur{i}.x = xTe(:,:,i);
%             dTe_cur_cur{i} = testSVMatLevel(dTe_cur_cur{i},dTr_cur_cur{i});
%             dTe_cur_cur{i}.info = '@vox';
%             dTe_cur_cur{i}.infoDim = combFields;
%             newxTe = cat(2,newxTe,dTe_cur_cur{i}.yhat);
%         end
%         %Store data
%         dTr_cur.x = newxTr;
%         dTr_cur.info = curInfo;
%         dTr_cur.infoDim = 'vox';
%         dTr_cur.subLevel = dTr_cur_cur;
%         dTe_cur.x = newxTe;
%         dTe_cur.info = curInfo;
%         dTe_cur.infoDim = 'vox';
%         dTe_cur.subLevel = dTe_cur_cur;
%
%     case 'vox'
%         %Reorganize data
%         dTr_cur = dTr{1}; dTr_cur.x = [];
%         dTe_cur = dTe{1}; dTe_cur.x = [];
%         xTr = [];
%         xTe = [];
%         for i = 1:length(dTr)
%             xTr = cat(3,xTr,dTr{i}.x);
%             xTe = cat(3,xTe,dTe{i}.x);
%         end
%         newxTr = [];
%         newxTe = [];
%         for i = 1:size(xTr,3)
%             %Train
%             dTr_cur_cur{i} = dTr{i};
%             dTr_cur_cur{i}.x = xTr(:,:,i);
%             dTr_cur_cur{i} = trainSVMatLevel2(dTr_cur_cur{i},p);
%             newxTr = cat(2,newxTr,dTr_cur_cur{i}.yhat);
%             %Test
%             dTe_cur_cur{i} = dTe{i};
%             dTe_cur_cur{i}.x = xTe(:,:,i);
%             dTe_cur_cur{i} = testSVMatLevel2(dTe_cur_cur{i},dTr_cur_cur{i});
%             newxTe = cat(2,newxTe,dTe_cur_cur{i}.yhat);
%         end
%         %Store data
%         dTr_cur.x = newxTr;
%         dTr_cur.info = curInfo;
%         dTr_cur.infoDim = combFields;
%         dTr_cur.subLevel = dTr_cur_cur;
%         dTe_cur.x = newxTe;
%         dTe_cur.info = curInfo;
%         dTe_cur.infoDim = combFields;
%         dTe_cur.subLevel = dTe_cur_cur;
%
%     case {'infoOrVox','voxOrInfo'}
%         error('double-check that')
%         %Reorganize data
%         dTr_cur = dTr{1}; dTr_cur.x = [];
%         dTe_cur = dTe{1}; dTe_cur.x = [];
%         newxTr = [];
%         newxTe = [];
%         for i = 1:length(dTr)
%             %Train
%             dTr{i} = trainSVMatLevel(dTr{i},p);
%             newxTr = cat(2,newxTr,dTr{i}.yhat);
%             %Test
%             dTe{i} = testSVMatLevel(dTe{i},dTr{i});
%             newxTe = cat(2,newxTe,dTe{i}.yhat);
%         end
%         %Store data
%         dTr_cur.x = newxTr;
%         dTr_cur.info = curInfo;
%         dTr_cur.infoDim = combFields;
%         dTr_cur.subLevel = dTr;
%         dTe_cur.x = newxTe;
%         dTe_cur.info = curInfo;
%         dTe_cur.infoDim = combFields;
%         dTe_cur.subLevel = dTe;
%
%     otherwise
%         error('X')
% end
%
%

function d = trainSVMatLevel2(d,p)
%% At fold
if ~p.subCrossVal.train || ~p.subCrossVal.test
    %Normalize
    y = d.y;
    x = d.x;
    norm.shift = mean(x,1);
    x = x-repmat(norm.shift,size(x,1),1);
    norm.scale = std(x,[],1);
    x = x./repmat(norm.scale,size(x,1),1);
    %Train model
    svmStruct = svmtrain(y,x,'-s 0 -t 0 -c 1 -q');
    if ~p.subCrossVal.test
        d.norm = norm;
        d.svmStruct = svmStruct;
    end
    if ~p.subCrossVal.train
        %Apply model to each data point
        [pred,acc,yhat] = svmpredict(y,x,svmStruct,'-q');
        %Extract model weigths
        [w,a] = getAandW(svmStruct,x);
        %Put in struct
        d.pred = pred;
        d.yhat = yhat;
        d.acc = acc(1);
        d.w = w;
        d.a = a;
    end
end

%% At subFold
if p.subCrossVal.train || p.subCrossVal.test
    crossVal = unique(d.crossVal);
    for kInd = 1:length(crossVal)
        
        %Train
        trInd = d.crossVal~=crossVal(kInd);
        dTr = [];
        dTr.x = d.x(trInd,:);
        dTr.y = d.y(trInd,:);
        pTr = p;
        pTr.subCrossVal.train = 0;
        pTr.subCrossVal.test = 0;
        dTr = trainSVMatLevel2(dTr,pTr);
        if p.subCrossVal.test
            d.norm(kInd,1) = dTr.norm;
            d.svmStruct(kInd,1) = dTr.svmStruct;
        end
        if p.subCrossVal.train
            d.w(kInd,:) = dTr.w;
            d.a(kInd,:) = dTr.a;
        end
        
        %Test
        if p.subCrossVal.train
            teInd = d.crossVal==crossVal(kInd);
            dTe = [];
            dTe.x = d.x(teInd,:);
            dTe.y = d.y(teInd,:);
            dTe = testSVMatLevel2(dTe,dTr);
            
            d.yhat(d.crossVal == crossVal(kInd),1) = dTe.yhat;
            d.pred(d.crossVal == crossVal(kInd),1) = dTe.pred;
            d.acc(kInd,1) = dTe.acc;
        end
    end
    d.acc = mean(d.acc,1);
    d.w = mean(d.w,1);
    d.a = mean(d.a,1);
end


function dTe = testSVMatLevel2(dTe,dTr)
if length(dTr.svmStruct)>1
    for subFoldInd = 1:length(dTr.svmStruct)
        %Normalize (at each subFold)
        x = dTe.x-repmat(dTr.norm(subFoldInd,1).shift,size(dTe.x,1),1);
        x = x./repmat(dTr.norm(subFoldInd,1).scale,size(x,1),1);
        %Test model (at each subFold)
        [dTe.pred(:,subFoldInd),acc,dTe.yhat(:,subFoldInd)] = svmpredict(dTe.y,x,dTr.svmStruct(subFoldInd,1),'-q');
        dTe.acc = acc(1);
    end
%     %Average across subFolds
%     dTe.yhat = mean(dTe.yhat,2);
%     %and recompute predictions and accuracies
%     dTe.pred = (dTe.yhat<0)+1;
%     dTe.acc = length(find(dTe.pred==dTe.y))/length(dTe.y)*100;
else
    %Normalize
    x = dTe.x-repmat(dTr.norm.shift,size(dTe.x,1),1);
    x = x./repmat(dTr.norm.scale,size(x,1),1);
    %Test model
    [dTe.pred(:,1),acc,dTe.yhat(:,1)] = svmpredict(dTe.y,x,dTr.svmStruct,'-q');
    dTe.acc = acc(1);
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


