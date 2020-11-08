function [rTr,rTe] = runMultiLevelSVM2(dOrig,p)
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
for kInd = 1:length(unique(d.crossVal))
    tic
    display(['Fold ' num2str(kInd) '/' num2str(length(unique(d.crossVal)))])
    
    %% Split train and test trials
    teInd = find(d.crossVal==kInd);
    trInd = find(d.crossVal~=kInd);
    
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
        case {'ALrev','ALrev_cross1','ALrev_cross2'}
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
            for infoInd = 1:length(info)
                curInfo = [info{infoInd} '__' p.infoComb];
                combInfo = info(infoInd);
                combInfoInd = ismember(info,combInfo);
                [dTr2{end+1}, dTe2{end+1}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),p,combInfo,curInfo,'vox');
                info{end+1} = curInfo;
            end
            
            
            
            
            
            % Collapse all info
            curInfo = p.infoComb;
            combInfo = info;
            combInfoInd = ismember(info,combInfo);
            [dTr3{kInd}, dTe3{kInd}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),p,combInfo,curInfo,'infoOrVox');
            dTr3{kInd} = trainSVMatLevel(dTr3{kInd},p);
            dTe3{kInd} = testSVMatLevel(dTe3{kInd},dTr3{kInd});
            clear dTr2 dTe2
            
        case {'all','all_cross1','all_cross2'}
            %Store voxel vector for each information type
            info = fields(tr);
            for infoInd = 1:length(info)
                dTr2{infoInd} = tr.(info{infoInd});
                dTr2{infoInd}.info = info{infoInd};
                dTr2{infoInd}.infoDim = 'vox';
                
                dTe2{infoInd} = te.(info{infoInd});
                dTe2{infoInd}.info = info{infoInd};
                dTe2{infoInd}.infoDim = 'vox';
            end
            clear tr te
            
            %Collapse information at the voxel level
            curInfo = 'delayCart__infoThenVox';
            combInfo = {'delayCartr' 'delayCarti'};
            combInfoInd = ismember(info,combInfo);
            [dTr2{end+1}, dTe2{end+1}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),p,combInfo,curInfo,'info');
            info{end+1} = curInfo;
            
            curInfo = 'delayCart_roiDelay__infoThenVox';
            combInfo = {'delayCartr_roiDelay' 'delayCarti_roiDelay'};
            combInfoInd = ismember(info,combInfo);
            [dTr2{end+1}, dTe2{end+1}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),p,combInfo,curInfo,'info');
            info{end+1} = curInfo;
            
            curInfo = 'cart_roiDelay__infoThenVox';
            combInfo = {'cartr_roiDelay' 'carti_roiDelay'};
            combInfoInd = ismember(info,combInfo);
            [dTr2{end+1}, dTe2{end+1}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),p,combInfo,curInfo,'info');
            info{end+1} = curInfo;
            
            curInfo = 'pol__infoThenVox';
            combInfo = {'delayPol' 'magPol'};
            combInfoInd = ismember(info,combInfo);
            [dTr2{end+1}, dTe2{end+1}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),p,combInfo,curInfo,'info');
            info{end+1} = curInfo;
                        
            
            % Collapse some information together
            curInfo = 'delayCart__voxThenInfo';
            combInfo = {'delayCartr' 'delayCarti'};
            combInfoInd = ismember(info,combInfo);
            [dTr2{end+1}, dTe2{end+1}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),p,combInfo,curInfo,'vox');
            info{end+1} = curInfo;
            
            curInfo = 'delayCart_roiDelay__voxThenInfo';
            combInfo = {'delayCartr_roiDelay' 'delayCarti_roiDelay'};
            combInfoInd = ismember(info,combInfo);
            [dTr2{end+1}, dTe2{end+1}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),p,combInfo,curInfo,'vox');
            info{end+1} = curInfo;
            
            curInfo = 'cart_roiDelay__voxThenInfo';
            combInfo = {'cartr_roiDelay' 'carti_roiDelay'};
            combInfoInd = ismember(info,combInfo);
            [dTr2{end+1}, dTe2{end+1}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),p,combInfo,curInfo,'vox');
            info{end+1} = curInfo;
            
            curInfo = 'pol__voxThenInfo';
            combInfo = {'delayPol' 'magPol'};
            combInfoInd = ismember(info,combInfo);
            [dTr2{end+1}, dTe2{end+1}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),p,combInfo,curInfo,'vox');
            info{end+1} = curInfo;
            
            
            % Collapse all info
            curInfo = p.infoComb;
            combInfo = info;
            combInfoInd = ismember(info,combInfo);
            [dTr3{kInd}, dTe3{kInd}] = svmReduction(dTr2(combInfoInd),dTe2(combInfoInd),p,combInfo,curInfo,'infoOrVox');
            dTr3{kInd} = trainSVMatLevel(dTr3{kInd},p);
            dTe3{kInd} = testSVMatLevel(dTe3{kInd},dTr3{kInd});
            clear dTr2 dTe2
        otherwise
            error('X')
    end
    toc
    
                
end


%% Remove unnecessary data
for kInd = 1:length(unique(d.crossVal))
    dTr3{kInd} = removex(dTr3{kInd});
    dTe3{kInd} = removex(dTe3{kInd});
end
%% Compile folds
rTr = compileFolds(dTr3,d);
rTe = compileFolds(dTe3,d);





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
        r.norm(k) = dFold{k}.norm;
    end
    if isfield(r,'svmStruct')
        r.svmStruct(:,k) = dFold{k}.svmStruct;
    end
    ind = ismember(d.crossVal,dFold{k}.crossVal);
    r.pred(ind,k) = dFold{k}.pred;
    r.yhat(ind,k) = dFold{k}.yhat;
    if isfield(r,'acc')
        r.acc(1,k) = dFold{k}.acc*length(dFold{k}.y);
    end
    if isfield(r,'w')
        r.w(k,:) = dFold{k}.w;
        r.a(k,:) = dFold{k}.a;
    end
end
if isfield(r,'acc')
    r.acc = sum(r.acc)./length(r.y);
%     r.acc = r.acc./(length(r.y)/length(dFold));
end
if isfield(r,'w')
    r.w= mean(r.w,1);
    r.a= mean(r.a,1);
end
r.yhat = nanmean(r.yhat,2);
r.pred = nanmean(r.pred,2);

if isfield(r,'subLevel')
    for subInd = 1:length(dFold{1}.subLevel)
        dFold2 = cell(1,length(dFold));
        for k = 1:length(dFold)
            dFold2{k} = dFold{k}.subLevel{subInd};
        end
        r.subLevel{subInd} = compileFolds(dFold2,d);
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
            dTr_cur_cur{i} = trainSVMatLevel(dTr_cur_cur{i},p);
            dTr_cur_cur{i}.info = '@vox';
            dTr_cur_cur{i}.infoDim = combFields;
            newxTr = cat(2,newxTr,dTr_cur_cur{i}.yhat);
            %Test
            dTe_cur_cur{i} = dTe_cur;
            dTe_cur_cur{i}.x = xTe(:,:,i);
            dTe_cur_cur{i} = testSVMatLevel(dTe_cur_cur{i},dTr_cur_cur{i});
            dTe_cur_cur{i}.info = '@vox';
            dTe_cur_cur{i}.infoDim = combFields;
            newxTe = cat(2,newxTe,dTe_cur_cur{i}.yhat);
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
            dTr_cur_cur{i} = trainSVMatLevel(dTr_cur_cur{i},p);
            newxTr = cat(2,newxTr,dTr_cur_cur{i}.yhat);
            %Test
            dTe_cur_cur{i} = dTe{i};
            dTe_cur_cur{i}.x = xTe(:,:,i);
            dTe_cur_cur{i} = testSVMatLevel(dTe_cur_cur{i},dTr_cur_cur{i});
            newxTe = cat(2,newxTe,dTe_cur_cur{i}.yhat);
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
            dTr{i} = trainSVMatLevel(dTr{i},p);
            newxTr = cat(2,newxTr,dTr{i}.yhat);
            %Test
            dTe{i} = testSVMatLevel(dTe{i},dTr{i});
            newxTe = cat(2,newxTe,dTe{i}.yhat);
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


function [dTr_cur,dTe_cur] = svmReduction2(dTr,dTe,p,combFields,curInfo,redDim)

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
            dTr_cur_cur{i} = trainSVMatLevel(dTr_cur_cur{i},p);
            dTr_cur_cur{i}.info = '@vox';
            dTr_cur_cur{i}.infoDim = combFields;
            newxTr = cat(2,newxTr,dTr_cur_cur{i}.yhat);
            %Test
            dTe_cur_cur{i} = dTe_cur;
            dTe_cur_cur{i}.x = xTe(:,:,i);
            dTe_cur_cur{i} = testSVMatLevel(dTe_cur_cur{i},dTr_cur_cur{i});
            dTe_cur_cur{i}.info = '@vox';
            dTe_cur_cur{i}.infoDim = combFields;
            newxTe = cat(2,newxTe,dTe_cur_cur{i}.yhat);
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
            newxTe = cat(2,newxTe,dTe_cur_cur{i}.yhat);
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
        error('double-check that')
        %Reorganize data
        dTr_cur = dTr{1}; dTr_cur.x = [];
        dTe_cur = dTe{1}; dTe_cur.x = [];
        newxTr = [];
        newxTe = [];
        for i = 1:length(dTr)
            %Train
            dTr{i} = trainSVMatLevel(dTr{i},p);
            newxTr = cat(2,newxTr,dTr{i}.yhat);
            %Test
            dTe{i} = testSVMatLevel(dTe{i},dTr{i});
            newxTe = cat(2,newxTe,dTe{i}.yhat);
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



function d = trainSVMatLevel(d,p,skipZforDelayCart)
if ~exist('skipZforDelayCart','var'); skipZforDelayCart = 0; end
%% Initiate stuff
d.yhat = nan(size(d.y));
d.pred = nan(size(d.y));
crossVal = unique(d.crossVal);
if p.getPattern
    if p.subCrossVal
        d.a = nan(length(crossVal),size(d.x,2));
        d.w = nan(length(crossVal),size(d.x,2));
    else
        d.a = nan(1,size(d.x,2));
        d.w = nan(1,size(d.x,2));
    end
end

%% Do on subFolds
if p.subCrossVal
    %Redefine crossVal
    d.crossVal_orig = d.crossVal;
    for kInd = 1:length(crossVal)
        d.crossVal(d.crossVal==crossVal(kInd)) = kInd;
    end
    clear crossVal
    
    %Loop over cross-validation folds
    crossVal = unique(d.crossVal);
    for kInd = 1:length(crossVal)
        trInd = d.crossVal~=kInd;
        teInd = d.crossVal==kInd;
        dTr = d;
        dTr.x = d.x(trInd,:);
        dTr.y = d.y(trInd,:);
        dTr.crossVal = d.crossVal(trInd,:);
        dTe = d;
        dTe.x = d.x(teInd,:);
        dTe.y = d.y(teInd,:);
        dTe.crossVal = d.crossVal(teInd,:);
        
        %Normalize
        if ~skipZforDelayCart
            dTr.norm.shift = mean(dTr.x,1);
        else
            dTr.norm.shift = zeros(1,size(dTr.x,2));
        end
        dTr.x = dTr.x-repmat(dTr.norm.shift,size(dTr.x,1),1);
        dTe.x = dTe.x-repmat(dTr.norm.shift,size(dTe.x,1),1);
        if ~skipZforDelayCart
            dTr.norm.scale = std(dTr.x,[],1);
        else
            dTr.norm.scale = ones(1,size(dTr.x,2));
        end
        dTr.x = dTr.x./repmat(dTr.norm.scale,size(dTr.x,1),1);
        dTe.x = dTe.x./repmat(dTr.norm.scale,size(dTe.x,1),1);
        %Train model
        dTr.svmStruct = svmtrain(dTr.y,dTr.x,['-s 0 -t 0 -c 1 -q']);
        %Test model
        [d.pred(teInd), ~, d.yhat(teInd)] = svmpredict(dTe.y,dTe.x,dTr.svmStruct,'-q');
        %Extract model weigths
        if p.getPattern
            [d.w(kInd,:),d.a(kInd,:)] = getAandW(dTr.svmStruct,dTr.x);
        end
        %Store subLevels into level
        if p.subCrossVal==2
            d.svmStruct(kInd,1) = dTr.svmStruct;
            d.norm.shift(kInd,:) = dTr.norm.shift;
            d.norm.scale(kInd,:) = dTr.norm.scale;
        end
    end
    if p.getPattern
        d.w = mean(d.w,1);
        d.a = mean(d.a,1);
    end
    d.crossVal = d.crossVal_orig; d = rmfield(d,'crossVal_orig');
    %Store subLevels into level
%     if p.subCrossVal==2
%         d.norm.shift = mean(d.norm.shift,1);
%         d.norm.scale = mean(d.norm.scale,1);
%     end
end

%% Do on all data
if p.subCrossVal==0 || p.subCrossVal==1
    %Normalize
    x = d.x;
    if ~skipZforDelayCart
        d.norm.shift = mean(x,1);
    else
        d.norm.shift = zeros(1,size(x,2));
    end
    x = x-repmat(d.norm.shift,size(x,1),1);
    if ~skipZforDelayCart
        d.norm.scale = std(x,[],1);
    else
        d.norm.scale = ones(1,size(x,2));
    end
    x = x./repmat(d.norm.scale,size(x,1),1);
    %Train model
    d.svmStruct = svmtrain(d.y,x,['-s 0 -t 0 -c 1 -q']);
    if p.subCrossVal==0
        %Apply model to each data point
        [d.pred, ~, d.yhat] = svmpredict(d.y,x,d.svmStruct,'-q');
        %Extract model weigths
        if p.getPattern
            [d.w,d.a] = getAandW(d.svmStruct,x);
        end
    end
end




% function d = trainSVMatLevel2(d,p,skipZforDelayCart)
% if ~exist('skipZforDelayCart','var'); skipZforDelayCart = 0; end
% %% Initiate stuff
% d.yhat = nan(size(d.y));
% d.pred = nan(size(d.y));
% crossVal = unique(d.crossVal);
% if p.getPattern
%     if p.subCrossVal
%         d.a = nan(length(crossVal),size(d.x,2));
%         d.w = nan(length(crossVal),size(d.x,2));
%     else
%         d.a = nan(1,size(d.x,2));
%         d.w = nan(1,size(d.x,2));
%     end
% end
% 
% %% Do on subFolds
% if p.subCrossVal
%     %Redefine crossVal
%     d.crossVal_orig = d.crossVal;
%     for kInd = 1:length(crossVal)
%         d.crossVal(d.crossVal==crossVal(kInd)) = kInd;
%     end
%     clear crossVal
%     
%     %Loop over cross-validation folds
%     crossVal = unique(d.crossVal);
%     for kInd = 1:length(crossVal)
%         trInd = d.crossVal~=kInd;
%         teInd = d.crossVal==kInd;
%         dTr = d;
%         dTr.x = d.x(trInd,:);
%         dTr.y = d.y(trInd,:);
%         dTr.crossVal = d.crossVal(trInd,:);
%         dTe = d;
%         dTe.x = d.x(teInd,:);
%         dTe.y = d.y(teInd,:);
%         dTe.crossVal = d.crossVal(teInd,:);
%         
%         %Normalize
%         if ~skipZforDelayCart
%             dTr.norm.shift = mean(dTr.x,1);
%         else
%             dTr.norm.shift = zeros(1,size(dTr.x,2));
%         end
%         dTr.x = dTr.x-repmat(dTr.norm.shift,size(dTr.x,1),1);
%         dTe.x = dTe.x-repmat(dTr.norm.shift,size(dTe.x,1),1);
%         if ~skipZforDelayCart
%             dTr.norm.scale = std(dTr.x,[],1);
%         else
%             dTr.norm.scale = ones(1,size(dTr.x,2));
%         end
%         dTr.x = dTr.x./repmat(dTr.norm.scale,size(dTr.x,1),1);
%         dTe.x = dTe.x./repmat(dTr.norm.scale,size(dTe.x,1),1);
%         %Train model
%         dTr.svmStruct = svmtrain(dTr.y,dTr.x,['-s 0 -t 0 -c 1 -q']);
%         %Test model
%         [d.pred(teInd), ~, d.yhat(teInd)] = svmpredict(dTe.y,dTe.x,dTr.svmStruct,'-q');
%         %Extract model weigths
%         if p.getPattern
%             [d.w(kInd,:),d.a(kInd,:)] = getAandW(dTr.svmStruct,dTr.x);
%         end
%         %Store subLevels into level
%         d.svmStruct(kInd,1) = dTr.svmStruct;
%         d.norm.shift(kInd,:) = dTr.norm.shift;
%         d.norm.scale(kInd,:) = dTr.norm.scale;        
%     end
%     d.w = mean(d.w,1);
%     d.a = mean(d.a,1);
%     d.norm.shift = mean(d.norm.shift,1);
%     d.norm.scale = mean(d.norm.scale,1);
%     d.crossVal = d.crossVal_orig; d = rmfield(d,'crossVal_orig');
% end





function dTe = testSVMatLevel(dTe,dTr)
if length(dTr.svmStruct)==1
    %Normalize
    dTe.x = dTe.x-repmat(dTr.norm.shift,size(dTe.x,1),1);
    dTe.x = dTe.x./repmat(dTr.norm.scale,size(dTe.x,1),1);
    %Test model
    [dTe.pred,acc,dTe.yhat] = svmpredict(dTe.y,dTe.x,dTr.svmStruct,'-q');
    dTe.acc = acc(1);
else
    for i = 1:size(dTr.svmStruct,1)
        %at every subFold models
        dTe.x = dTe.x-repmat(dTr.norm.shift(i,:),size(dTe.x,1),1);
        dTe.x = dTe.x./repmat(dTr.norm.scale(i,:),size(dTe.x,1),1);
        
        [dTe.pred(:,i),~,dTe.yhat(:,i)] = svmpredict(dTe.y,dTe.x,dTr.svmStruct(i,1),'-q');
    end
    %the average
    dTe.yhat = mean(dTe.yhat,2);
    %and recompute predictions and accuracies
    dTe.pred = (dTe.yhat<0)+1;
    dTe.acc = length(find(dTe.pred==dTe.y))/length(dTe.y)*100;
end

function dTe = testSVMatLevel2(dTe,dTr)
%Normalize
dTe.x = dTe.x-repmat(dTr.norm.shift,size(dTe.x,1),1);
dTe.x = dTe.x./repmat(dTr.norm.scale,size(dTe.x,1),1);
%Test model
for i = 1:size(dTr.svmStruct,1)
    %at every subFold models
    [dTe.pred(:,i),~,dTe.yhat(:,i)] = svmpredict(dTe.y,dTe.x,dTr.svmStruct(i,1),'-q');
end
%the average
dTe.yhat = mean(dTe.yhat,2);
%and recompute predictions and accuracies
dTe.pred = (dTe.yhat<0)+1;
dTe.acc = length(find(dTe.pred==dTe.y))/length(dTe.y)*100;


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



% function dTr_cur = trainAtVox(dTr,p,combFields,curInfo)
% dTr_cur = dTr{1}; dTr_cur.x = [];
% for i = 1:length(dTr)
%     dTr_cur.x = cat(3,dTr_cur.x,dTr{i}.x);
% end
% dTr_cur.x = permute(dTr_cur.x,[1 3 2]);
% newx = [];
% for vox = 1:size(dTr_cur.x,3)
%     dTr_cur_cur{vox} = dTr_cur;
%     dTr_cur_cur{vox}.x = dTr_cur.x(:,:,vox);
%     dTr_cur_cur{vox} = trainSVMatLevel(dTr_cur_cur{vox},p);
%     dTr_cur_cur{vox}.info = combFields;
%     %                 dTr_cur_cur{vox}.info = ['vox' num2str(vox)];
%     newx = cat(2,newx,dTr_cur_cur{vox}.yhat);
% end
% dTr_cur.info = curInfo;
% dTr_cur.subLevel = dTr_cur_cur;
% dTr_cur.x = newx;


% function [d2Tr,d1Dim1Tr_ord1,d1Dim2Tr_ord1,d1Dim1Tr_ord2,d1Dim2Tr_ord2,d2Te,d1Dim1Te_ord1,d1Dim2Te_ord1,d1Dim1Te_ord2,d1Dim2Te_ord2,curInd,curInd1,curInd2] = trainNtest_2levels(d,p,dTr,dTe,skipZforDelayCart)
% if ~exist('skipZforDelayCart','var')
%     skipZforDelayCart = 0;
% end
% %% Loop over label pairs
% %%%%%
% %%%%%
% %% Train
% %%%%%
% %%%%%
% %     dDim1Tr_ord1 = cell(size(d.labelPairs,1),sz(end));% dDim1Te = cell(size(d1Tr.x,3),1);
% %     dDim2Tr = cell(size(d.labelPairs,1),1);
% for pairInd = 1:size(d.labelPairs,1)
%     
%     %% Hierarchical
%     % Prepare data
%     [d1Tr,~,~,~] = getLabelPairData(dTr,d.labelPairs(pairInd,:));
%     
%     % Run SVM on freq then chan
%     [d1Dim1Tr_ord1(pairInd,:),d1Dim2Tr_ord1(pairInd)] = train2Levels(d1Tr,p);
%     
%     % Run SVM on chan then freq
%     d1Tr.x = permute(d1Tr.x,[1 3 2]);
%     %         d1Tr.x = cat(1,d1Tr.x(end/2+1:end,:,:),d1Tr.x(1:end/2,:,:));
%     [d1Dim1Tr_ord2(pairInd,:),d1Dim2Tr_ord2(pairInd)] = train2Levels(d1Tr,p,[],skipZforDelayCart);
%     %         figure('WindowStyle','docked');
%     %         hist(d1Dim2Tr_ord2{1}.a); xlabel('a')
%     %         figure('WindowStyle','docked');
%     %         hist(d1Dim2Tr_ord2{1}.w); xlabel('w')
%     
%     %% Supra-Hierarchical
%     d2Tr(pairInd) = {d1Tr};
%     d2Tr{pairInd}.x = cat(2,d1Dim2Tr_ord1{pairInd}.yhat,d1Dim2Tr_ord2{pairInd}.yhat);
%     d2Tr{pairInd} = trainSVMatLevel(d2Tr{pairInd},p);
%     
%     %         %% Compile a and w
%     %         if p.getPattern
%     %             %From high to low levels
%     %             level3{1}.a(pairInd,:,kInd) = d2Tr{pairInd}.a;
%     %             level3{1}.w(pairInd,:,kInd) = d2Tr{pairInd}.w;
%     %             level3{1}.info = 'freqThenChan; chanThenFreq';
%     %             level2{1}.a(pairInd,:,kInd) = d1Dim2Tr_ord1{pairInd}.a;
%     %             level2{1}.w(pairInd,:,kInd) = d1Dim2Tr_ord1{pairInd}.w;
%     %             level2{1}.info(pairInd,:,kInd) = 'freqThenChan -> chan';
%     %             level2{2}.a(pairInd,:,kInd) = d1Dim2Tr_ord2{pairInd}.a;
%     %             level2{2}.w(pairInd,:,kInd) = d1Dim2Tr_ord2{pairInd}.w;
%     %             level2{2}.info(pairInd,:,kInd) = 'chanThenFreq -> freq';
%     %         end
% end
% 
% %%%%%
% %%%%%
% %% Test
% %%%%%
% %%%%%
% %     display('Testing all pairs')
% %     dDim1Te = cell(size(d.labelPairs,1),sz(end));% dDim1Te = cell(size(d1Tr.x,3),1);
% %     dDim2Te = cell(size(d.labelPairs,1),1);
% curInd1 = cell(1,size(d.labelPairs,1));
% curInd2 = cell(1,size(d.labelPairs,1));
% for pairInd = 1:size(d.labelPairs,1)
%     % Prepare data
%     [d1Te,curInd{pairInd},curInd1{pairInd},curInd2{pairInd}] = getLabelPairData(dTe,d.labelPairs(pairInd,:));
%     
%     %% Hierarchical
%     % SVM on Freq then Chanels
%     [d1Dim1Te_ord1(pairInd,:),d1Dim2Te_ord1(pairInd)] = test2Levels(d1Te,d1Dim1Tr_ord1(pairInd,:),d1Dim2Tr_ord1);
%     
%     % SVM on Channels then Freq
%     d1Te.x = permute(d1Te.x,[1 3 2]);
%     [d1Dim1Te_ord2(pairInd,:),d1Dim2Te_ord2(pairInd)] = test2Levels(d1Te,d1Dim1Tr_ord2(pairInd,:),d1Dim2Tr_ord2);
%     
%     
%     %% Supra-Hierarchical
%     d2Te(pairInd) = {d1Te};
%     d2Te{pairInd}.x = cat(2,d1Dim2Te_ord1{pairInd}.yhat,d1Dim2Te_ord2{pairInd}.yhat);
%     d2Te{pairInd} = testSVMatLevel(d2Te{pairInd},d2Tr{pairInd});
% end

% function [tr,te,curInd,curInd1,curInd2] = trainNtest_L1L2(d,p,dTr,dTe,signOutput)
% if ~exist('signOutput','var')
%     signOutput = 0;
% end
% %% Loop over label pairs
% %%%%%
% %% Train
% %%%%%
% tr.L1 = cell(size(d.labelPairs,1),size(dTr.x,3));
% tr.L2 = cell(size(d.labelPairs,1),1);
% tr = repmat(tr,size(d.labelPairs,1),1);
% for pairInd = 1:size(d.labelPairs,1)
%     % Prepare data
%     [dTr_cur,~,~,~] = getLabelPairData(dTr,d.labelPairs(pairInd,:));
% %     curInd1 = dTr.y == d.labelPairs(pairInd,1); curInd2 = dTr.y == d.labelPairs(pairInd,2);
% %     curInd = curInd1|curInd2;
% %     dTr_cur = dTr;
% %     dTr_cur.y = cat(1,ones(length(find(curInd1)),1),ones(length(find(curInd2)),1).*2);
% %     dTr_cur.x = dTr.x(curInd,:,:);
% %     dTr_cur.crossVal = dTr.crossVal(curInd);
% %     dTr_cur.labelPair = d.labelPairs(pairInd,:);
%     % Train two levels
%     [trL1,trL2] = train2Levels(dTr_cur,p,signOutput);
%     tr(pairInd).L1 = trL1;
%     tr(pairInd).L2 = trL2;
%     
%     %         %% Compile a and w
%     %         if p.getPattern
%     %             %From high to low levels
%     %             level3{1}.a(pairInd,:,kInd) = d2Tr{pairInd}.a;
%     %             level3{1}.w(pairInd,:,kInd) = d2Tr{pairInd}.w;
%     %             level3{1}.info = 'freqThenChan; chanThenFreq';
%     %             level2{1}.a(pairInd,:,kInd) = d1Dim2Tr_ord1{pairInd}.a;
%     %             level2{1}.w(pairInd,:,kInd) = d1Dim2Tr_ord1{pairInd}.w;
%     %             level2{1}.info(pairInd,:,kInd) = 'freqThenChan -> chan';
%     %             level2{2}.a(pairInd,:,kInd) = d1Dim2Tr_ord2{pairInd}.a;
%     %             level2{2}.w(pairInd,:,kInd) = d1Dim2Tr_ord2{pairInd}.w;
%     %             level2{2}.info(pairInd,:,kInd) = 'chanThenFreq -> freq';
%     %         end
% end
% 
% %%%%%
% %% Test
% %%%%%
% te.L1 = cell(size(d.labelPairs,1),size(dTr.x,3));
% te.L2 = cell(size(d.labelPairs,1),1);
% te = repmat(te,size(d.labelPairs,1),1);
% curInd1 = cell(1,size(d.labelPairs,1));
% curInd2 = cell(1,size(d.labelPairs,1));
% for pairInd = 1:size(d.labelPairs,1)
%     % Prepare data
%     [dTe_cur,curInd{pairInd},curInd1{pairInd},curInd2{pairInd}] = getLabelPairData(dTe,d.labelPairs(pairInd,:));
% %     curInd1{pairInd} = dTe.y == d.labelPairs(pairInd,1); curInd2{pairInd} = dTe.y == d.labelPairs(pairInd,2);
% %     curInd = curInd1{pairInd}|curInd2{pairInd};
% %     dTe_cur = dTe;
% %     dTe_cur.y = cat(1,ones(length(find(curInd1{pairInd})),1),ones(length(find(curInd2{pairInd})),1).*2);
% %     dTe_cur.x = dTe.x(curInd,:,:);
% %     dTe_cur.labelPair = d.labelPairs(pairInd,:);
%     
%     % Test 2 levels
%     [teL1,teL2] = test2Levels(dTe_cur,trL1(pairInd,:),trL2(pairInd),signOutput);
%     te(pairInd).L1 = teL1;
%     te(pairInd).L2 = teL2;
% end


% function [dDim1Tr,dDim2Tr] = train2Levels(d1Tr,p,signOutput,skipZforDelayCart)
% if ~exist('signOutput','var')|isempty(signOutput); signOutput = 0; end
% if ~exist('skipZforDelayCart','var'); skipZforDelayCart = 0; end
% sz = size(d1Tr.x);
% %% First-level
% dDim1Tr = cell(1,sz(end));
% dDim2Tr = cell(1,1);
% for chanInd = 1:sz(end)
%     %Get data
%     dDim1Tr{chanInd} = d1Tr;
%     dDim1Tr{chanInd}.x = d1Tr.x(:,:,chanInd);
%     %Do SVM
%     if signOutput; %Sign the output acording to the sign of the imaginary value
%         pCur = p; pCur.getPattern = 1;
%         dDim1Tr{chanInd} = trainSVMatLevel(dDim1Tr{chanInd},pCur,skipZforDelayCart);
%         if dDim1Tr{chanInd}.w(2)~=0
%             dDim1Tr{chanInd}.yhat = sign(dDim1Tr{chanInd}.w(2)).*dDim1Tr{chanInd}.yhat;
%         end
%     else
%         dDim1Tr{chanInd} = trainSVMatLevel(dDim1Tr{chanInd},p,skipZforDelayCart);
%     end
% end
% 
% %% Second-level
% %Get data
% if dDim1Tr{1}.getPattern
%     dDim2Tr{1} = rmfield(dDim1Tr{1},{'norm' 'svmStruct' 'yhat' 'bias' 'w' 'a'});%cell(size(d.x,3),1); dDim1Te = cell(size(d.x,3),1);
% else
%     dDim2Tr{1} = rmfield(dDim1Tr{1},{'norm' 'svmStruct' 'yhat' 'bias'});%cell(size(d.x,3),1); dDim1Te = cell(size(d.x,3),1);
% end
% dDim2Tr{1}.x = nan(size(dDim1Tr{1}.yhat,1),sz(end));
% for chanInd = 1:length(dDim1Tr)
%     dDim2Tr{1}.x(:,chanInd) = dDim1Tr{chanInd}.yhat;
% end
% %Do SVM
% dDim2Tr{1} = trainSVMatLevel(dDim2Tr{1},p);

% function [dDim1Te,dDim2Te] = test2Levels(d1Te,dDim1Tr,dDim2Tr,signOutput)
% if ~exist('signOutput','var'); signOutput = 0; end
% sz = size(d1Te.x);
% dDim1Te = cell(1,sz(end));
% dDim2Te = cell(1,1);
% %% First-level
% for chanInd = 1:sz(end)
%     %Get data
%     dDim1Te{chanInd} = d1Te;
%     dDim1Te{chanInd}.x = d1Te.x(:,:,chanInd);
%     %Do SVM
%     dDim1Te{chanInd} = testSVMatLevel(dDim1Te{chanInd},dDim1Tr{chanInd});
%     if signOutput
%         if dDim1Tr{chanInd}.w(2)~=0
%             dDim1Te{chanInd}.yhat = dDim1Te{chanInd}.yhat.*sign(dDim1Tr{chanInd}.w(2));
%         end
%     end
% end
% 
% %% Second-level: SVM on chanel dimension, using output from first-level
% %Get data
% dDim2Te{1} = rmfield(dDim1Te{chanInd},{'yhat' 'acc'});
% dDim2Te{1}.x = nan(size(dDim1Te{chanInd}.yhat,1),sz(end));
% for chanInd = 1:length(dDim1Tr)
%     dDim2Te{1}.x(:,chanInd) = dDim1Te{chanInd}.yhat;
% end
% %Do SVM
% dDim2Te{1} = testSVMatLevel(dDim2Te{1},dDim2Tr{1});




% function [d,curInd,curInd1,curInd2] = getLabelPairData(d,labelPairs)
% curInd1 = d.y == labelPairs(1); curInd2 = d.y == labelPairs(2);
% curInd = curInd1|curInd2;
% d.y = cat(1,ones(length(find(curInd1)),1),ones(length(find(curInd2)),1).*2);
% d.x = d.x(curInd,:,:);
% d.crossVal = d.crossVal(curInd);
% d.labelPair = labelPairs;
