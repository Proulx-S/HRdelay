function [dTr4,dVal4,dTe4] = runMultiLevelSVM6(dOrig,p,verbose)

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
%     else
%         fprintf ('Fold %d/%d; \r',crossVal(kInd),length(crossVal));
    end
    switch p.infoComb
        case {'catComb_subFold' 'catComb_subFold_AL'}
            if strcmp(p.infoComb(end-1:end),'AL'); ALflag = 1; else ALflag = 0; end
            subCrossVal = crossVal(crossVal~=crossVal(kInd));
            for kSubInd = 1:length(subCrossVal)
                %% Split train and test trials
                [dTr{1,1},dVal{1,1},dTe{1,1},trInd{1,kSubInd},valInd{1,kSubInd},teInd{1,kSubInd},info] = splitTrainValTest(d,crossVal(kInd),subCrossVal(kSubInd));
                
                %% Concat info
                info = reshape(info,[2 length(info)/2])';
                infoComb = {'pol' 'delayCart' 'delayCart_roiDelay' 'cart' 'cart_roiDelay'}';
                for infoInd = 1:size(info,1)
                    [dTr_tmp,dVal_tmp,dTe_tmp] = catInfo(d,dTr(1,1),dVal(1,:),dTe(1,1),trInd(1,1),valInd(1,1),teInd(1,1),info(infoInd,:),infoComb{infoInd});
                    dTr2{1,kSubInd}(infoInd) = dTr_tmp{1,1};
                    dVal2{1,kSubInd}(infoInd) = dVal_tmp{1,1};
                    dTe2{1,kSubInd}(infoInd) = dTe_tmp{1,1};
                end
                clear dTr dVal dTe clear dTr_tmp dVal_tmp dTe_tmp
                [dTr2{1,kSubInd},dVal2{1,kSubInd},dTe2{1,kSubInd}] = trainNtest(dTr2{1,kSubInd},dVal2{1,kSubInd},dTe2{1,kSubInd},p,ALflag);
            end
            
            %% Bring up one level
            [dTr3(kInd,:),dVal3(kInd,:),dTe3(kInd,:)] = levelUp(d,dTr2(1,:),dVal2(1,:),dTe2(1,:),trInd(1,:),valInd(1,:),teInd(1,:),[],p.infoComb);
            
            
            %% SVM on higher level
            for kSubInd = 1:length(subCrossVal)
                [dTr3{kInd,kSubInd},dVal3{kInd,kSubInd},dTe3{kInd,kSubInd}] = trainNtest(dTr3{kInd,kSubInd},dVal3{kInd,kSubInd},dTe3{kInd,kSubInd},p,0);
            end
        case {'infoComb_subFold' 'infoComb_subFold_AL'}
            if strcmp(p.infoComb(end-1:end),'AL'); ALflag = 1; else ALflag = 0; end
            subCrossVal = crossVal(crossVal~=crossVal(kInd));
            for kSubInd = 1:length(subCrossVal)
                
                %% Split train and test trials
                [dTr{1,kSubInd},dVal{1,kSubInd},dTe{1,kSubInd},trInd{1,kSubInd},valInd{1,kSubInd},teInd{1,kSubInd},info] = splitTrainValTest(d,crossVal(kInd),subCrossVal(kSubInd));
                
                %% SVM on low-level (collapse voxels)
                [dTr{1,kSubInd},dVal{1,kSubInd},dTe{1,kSubInd}] = trainNtest(dTr{1,kSubInd},dVal{1,kSubInd},dTe{1,kSubInd},p,ALflag);
            end
            
            %% Collapse info
            info = reshape(info,[2 length(info)/2])';
            infoComb = {'pol' 'delayCart' 'delayCart_roiDelay' 'cart' 'cart_roiDelay'}';
            for infoInd = 1:size(info,1)
                [dTr_tmp,dVal_tmp,dTe_tmp] = levelUp(d,dTr(1,:),dVal(1,:),dTe(1,:),trInd(1,:),valInd(1,:),teInd(1,:),info(infoInd,:),infoComb{infoInd});
                for kSubInd = 1:length(subCrossVal)
                    dTr2{1,kSubInd}(infoInd) = dTr_tmp{1,kSubInd};
                    dVal2{1,kSubInd}(infoInd) = dVal_tmp{1,kSubInd};
                    dTe2{1,kSubInd}(infoInd) = dTe_tmp{1,kSubInd};
                end
            end
            for kSubInd = 1:length(subCrossVal)
                [dTr2{1,kSubInd},dVal2{1,kSubInd},dTe2{1,kSubInd}] = trainNtest(dTr2{1,kSubInd},dVal2{1,kSubInd},dTe2{1,kSubInd},p,0);
            end
            
            %% SVM on higher level
            [dTr3(kInd,:),dVal3(kInd,:),dTe3(kInd,:)] = levelUp(d,dTr2(1,:),dVal2(1,:),dTe2(1,:),trInd(1,:),valInd(1,:),teInd(1,:),[],p.infoComb);
            for kSubInd = 1:length(subCrossVal)
                [dTr3{kInd,kSubInd},dVal3{kInd,kSubInd},dTe3{kInd,kSubInd}] = trainNtest(dTr3{kInd,kSubInd},dVal3{kInd,kSubInd},dTe3{kInd,kSubInd},p,0);
            end
            
            
            
        case {'infoComb_CsubFold'}
            if strcmp(p.infoComb(end-1:end),'AL'); ALflag = 1; else ALflag = 0; end
            subCrossVal = crossVal(crossVal~=crossVal(kInd));
            for kSubInd = 1:length(subCrossVal)
                kSubInd
                %% Split train and test trials
                [dTr{1,kSubInd},dVal{1,kSubInd},dTe{1,kSubInd},trInd{1,kSubInd},valInd{1,kSubInd},teInd{1,kSubInd},info] = splitTrainValTest(d,crossVal(kInd),subCrossVal(kSubInd));
                
                for cInd = 1:length(p.cList)
                    curP = p;
                    curP.C = p.cList(cInd);
                    %% SVM on low-level (collapse voxels)
                    [dTr_tmp{1,kSubInd,cInd},dVal_tmp{1,kSubInd,cInd},dTe_tmp{1,kSubInd,cInd}] = trainNtest(dTr{1,kSubInd},dVal{1,kSubInd},dTe{1,kSubInd},curP,ALflag);
                end
            end
            for cInd = 1:length(p.cList)
                dTr_tmp2(1,cInd) = collapseSubFold(d,dTr_tmp(1,kSubInd,cInd));
                dVal_tmp2(1,cInd) = collapseSubFold(d,dVal_tmp(1,kSubInd,cInd));
                dTe_tmp2(1,cInd) = collapseSubFold(d,dTe_tmp(1,kSubInd,cInd));
            end
            dTr_tmp2 = collapseSubFold(d,dTr_tmp2);
            dVal_tmp2 = collapseSubFold(d,dVal_tmp2);
            dTe_tmp2 = collapseSubFold(d,dTe_tmp2);
            
            keyboard
            figure('WindowStyle','docked');
            for anaInd = 1:length(dTr_tmp2{1})
                subplot(1,length(dTr_tmp2{1}),anaInd)
                distT_tr = [];
                dist_tr = [];
                dist_val = [];
                for cInd = 1:size(dTr_tmp2{1}{anaInd}.yhat,2)
                    [~,~,~,STATS] = ttest2(dTr_tmp2{1}{anaInd}.yhat(dTr_tmp2{1}{anaInd}.y(:,cInd)==1,cInd),dTr_tmp2{1}{anaInd}.yhat(dTr_tmp2{1}{anaInd}.y(:,cInd)==2,cInd));
                    distT_tr(cInd) = STATS.tstat;
                    [~,~,~,STATS] = ttest2(dVal_tmp2{1}{anaInd}.yhat(dVal_tmp2{1}{anaInd}.y(:,cInd)==1,cInd),dVal_tmp2{1}{anaInd}.yhat(dVal_tmp2{1}{anaInd}.y(:,cInd)==2,cInd));
                    distT_val(cInd) = STATS.tstat;

%                     dist_tr(cInd) = mean(dTr_tmp2{1}{anaInd}.yhat(dTr_tmp2{1}{anaInd}.y(:,cInd)==1,cInd),1) - mean(dTr_tmp2{1}{anaInd}.yhat(dTr_tmp2{1}{anaInd}.y(:,cInd)==2,cInd),1);
%                     dist_val(cInd) = mean(dVal_tmp2{1}{anaInd}.yhat(dVal_tmp2{1}{anaInd}.y(:,cInd)==1,cInd),1) - mean(dVal_tmp2{1}{anaInd}.yhat(dVal_tmp2{1}{anaInd}.y(:,cInd)==2,cInd),1);
                end
                plot(log(p.cList),distT_tr,'-o'); hold on
                plot(log(p.cList),distT_val,'-o'); hold on

%                 plot(log(p.cList),dist_tr./max(dist_tr),'-o'); hold on
%                 plot(log(p.cList),dist_val./max(dist_val),'-o'); hold on
                title(dTr_tmp2{1}{anaInd}.info)
            end
            
            
            ax = gca;
ax.ColorOrderIndex = 1;
            
            
            
            
            figure('WindowStyle','docked');
            for anaInd = 1:length(dTr_tmp2{1})
                for cInd = 1:size(dTr_tmp2{1}{anaInd}.yhat,2)
                    plot(dTr_tmp2{1}{anaInd}.yhat(:,cInd)); hold on
                end
            end
            
            figure('WindowStyle','docked');
            for anaInd = 1:length(dTr_tmp2{1})
                plot(log(p.cList),dTr_tmp2{1}{anaInd}.acc); hold on
                plot(log(p.cList),dVal_tmp2{1}{anaInd}.acc,':'); hold on
            end
            
            %% Collapse info
            info = reshape(info,[2 length(info)/2])';
            infoComb = {'pol' 'delayCart' 'delayCart_roiDelay' 'cart' 'cart_roiDelay'}';
            for infoInd = 1:size(info,1)
                [dTr_tmp,dVal_tmp,dTe_tmp] = levelUp(d,dTr(1,:),dVal(1,:),dTe(1,:),trInd(1,:),valInd(1,:),teInd(1,:),info(infoInd,:),infoComb{infoInd});
                for kSubInd = 1:length(subCrossVal)
                    dTr2{1,kSubInd}(infoInd) = dTr_tmp{1,kSubInd};
                    dVal2{1,kSubInd}(infoInd) = dVal_tmp{1,kSubInd};
                    dTe2{1,kSubInd}(infoInd) = dTe_tmp{1,kSubInd};
                end
            end
            for kSubInd = 1:length(subCrossVal)
                [dTr2{1,kSubInd},dVal2{1,kSubInd},dTe2{1,kSubInd}] = trainNtest(dTr2{1,kSubInd},dVal2{1,kSubInd},dTe2{1,kSubInd},p,0);
            end
            
            %% SVM on higher level
            [dTr3(kInd,:),dVal3(kInd,:),dTe3(kInd,:)] = levelUp(d,dTr2(1,:),dVal2(1,:),dTe2(1,:),trInd(1,:),valInd(1,:),teInd(1,:),[],p.infoComb);
            for kSubInd = 1:length(subCrossVal)
                [dTr3{kInd,kSubInd},dVal3{kInd,kSubInd},dTe3{kInd,kSubInd}] = trainNtest(dTr3{kInd,kSubInd},dVal3{kInd,kSubInd},dTe3{kInd,kSubInd},p,0);
            end
            
        case {'voxComb_subFold' 'voxComb_subFold_AL'}
            if strcmp(p.infoComb(end-1:end),'AL'); ALflag = 1; else ALflag = 0; end
            subCrossVal = crossVal(crossVal~=crossVal(kInd));
            %% Split train and test trials... and collapsing some info on a voxel-by-voxel basis
            info = splitTrainValTest;
            info = reshape(info,[2 length(info)/2])';
            infoComb = {'pol' 'delayCart' 'delayCart_roiDelay' 'cart' 'cart_roiDelay'}';
            [dTr_tmp,dVal_tmp,dTe_tmp,trInd,valInd,teInd] = voxelSplitNcomb(d,p,crossVal(kInd),subCrossVal,info,infoComb);
            
            %% Collapse voxels at each info combination
            for kSubInd = 1:length(subCrossVal)                
                [dTr_tmp{1,kSubInd},dVal_tmp{1,kSubInd},dTe_tmp{1,kSubInd}] = trainNtest(dTr_tmp{1,kSubInd},dVal_tmp{1,kSubInd},dTe_tmp{1,kSubInd},p,ALflag);
            end
            [dTr_tmp,dVal_tmp,dTe_tmp] = levelUp(d,dTr_tmp,dVal_tmp,dTe_tmp,trInd,valInd,teInd,[],p.infoComb);
            
            %% Collapse all info
            for kSubInd = 1:length(subCrossVal)
                [dTr3{kInd,kSubInd},dVal3{kInd,kSubInd},dTe3{kInd,kSubInd}] = trainNtest(dTr_tmp{1,kSubInd},dVal_tmp{1,kSubInd},dTe_tmp{1,kSubInd},p,0);
            end

        case {'lowLevel_subFold' 'lowLevel_subFold_AL'}
            if strcmp(p.infoComb(end-1:end),'AL'); ALflag = 1; else ALflag = 0; end
            subCrossVal = crossVal(crossVal~=crossVal(kInd));
            for kSubInd = 1:length(subCrossVal)
                
                %% Split train and test trials
                [dTr{1,kSubInd},dVal{1,kSubInd},dTe{1,kSubInd},trInd{1,kSubInd},valInd{1,kSubInd},teInd{1,kSubInd},info] = splitTrainValTest(d,crossVal(kInd),subCrossVal(kSubInd));
                
                %% SVM on low-level
                [dTr2{1,kSubInd},dVal2{1,kSubInd},dTe2{1,kSubInd}] = trainNtest(dTr{1,kSubInd},dVal{1,kSubInd},dTe{1,kSubInd},p,ALflag);
                
            end
            
            %% Bring up one level
            [dTr3(kInd,:),dVal3(kInd,:),dTe3(kInd,:)] = levelUp(d,dTr2(1,:),dVal2(1,:),dTe2(1,:),trInd(1,:),valInd(1,:),teInd(1,:),info,p.infoComb);
            
            %% SVM on higher level
            for kSubInd = 1:length(subCrossVal)
                [dTr3{kInd,kSubInd},dVal3{kInd,kSubInd},dTe3{kInd,kSubInd}] = trainNtest(dTr3{kInd,kSubInd},dVal3{kInd,kSubInd},dTe3{kInd,kSubInd},p,0);
            end
        
        otherwise
            error('XX')
    end
end
% clear dTr dVal dTe dTr2 dVal2 dTe2
clearvars -except d dTr3 dVal3 dTe3

%% Rearange
for kInd = 1:size(dTr3,1)
    dTr4(kInd,1) = collapseSubFold(d,dTr3(kInd,:));
    dVal4(kInd,1) = collapseSubFold(d,dVal3(kInd,:));
    dTe4(kInd,1) = collapseSubFold(d,dTe3(kInd,:));
    dTr3(kInd,:) = cell(size(dTr3(kInd,:)));
    dVal3(kInd,:) = cell(size(dVal3(kInd,:)));
    dTe3(kInd,:) = cell(size(dTe3(kInd,:)));
end
clear dTr3 dVal3 dTe3

dTr4 = collapseFold(dTr4);
dVal4 = collapseFold(dVal4);
dTe4 = collapseFold(dTe4);

dTr4 = replaceD(d,dTr4,0); dTr4 = dTr4{1}{1};
dVal4 = replaceD(d,dVal4,0); dVal4 = dVal4{1}{1};
dTe4 = replaceD(d,dTe4); dTe4 = dTe4{1}{1};





function [dTr2,dVal2,dTe2,trInd,valInd,teInd,info] = splitTrainValTest(d,crossVal,subCrossVal)
if ~exist('d','var') && ~exist('crossVal','var') && ~exist('subCrossVal','var')
    dTr2 = {'magPol' 'delayPol' 'delayCartr' 'delayCarti' 'delayCartr_roiDelay' 'delayCarti_roiDelay' 'cartr' 'carti' 'cartr_roiDelay' 'carti_roiDelay'};
    dVal2 = [];
    dTe2 = [];
    trInd = [];
    valInd = [];
    teInd = [];
    info = [];
    return
end

trInd = find(d.crossVal~=crossVal & d.crossVal~=subCrossVal);
valInd = find(d.crossVal==subCrossVal);
teInd = find(d.crossVal==crossVal);

dTr = rmfield(d,{'label'});
dTr.x = d.x(trInd,:,:);
dTr.y = d.y(trInd,:,:);
dTr.crossVal = d.crossVal(trInd);
dVal = rmfield(d,{'label'});
dVal.x = d.x(valInd,:,:);
dVal.y = d.y(valInd,:,:);
dVal.crossVal = d.crossVal(valInd);
dTe = rmfield(d,{'label'});
dTe.x = d.x(teInd,:,:);
dTe.y = d.y(teInd,:,:);
dTe.crossVal = d.crossVal(teInd);

%% Normalize phase and mag
trPhase = angle(dTr.x);
tePhase = angle(dTe.x);
valPhase = angle(dVal.x);
trMag = abs(dTr.x);
teMag = abs(dTe.x);
valMag = abs(dVal.x);
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

val.magPol = dVal;
val.magPol.x = valMag./repmat(magScale,size(valMag,1),1);
val.delayPol = dVal;
val.delayPol.x = wrapToPi(valPhase-repmat(phaseShift,size(valPhase,1),1));
val.delayCartr = dVal; val.delayCarti = dVal;
[val.delayCartr.x,val.delayCarti.x] = pol2cart(wrapToPi(valPhase-repmat(phaseShift,size(valPhase,1),1)),ones(size(valMag)));
val.delayCartr_roiDelay = dVal; val.delayCarti_roiDelay = dVal;
[val.delayCartr_roiDelay.x,val.delayCarti_roiDelay.x] = pol2cart(wrapToPi(valPhase-repmat(phaseShift_roi,size(valPhase,1),1)),ones(size(valMag)));
val.cartr = dVal; val.carti = dVal;
[val.cartr.x,val.carti.x] = pol2cart(wrapToPi(valPhase-repmat(phaseShift,size(valPhase,1),1)),valMag./repmat(magScale,size(valMag,1),1));
val.cartr_roiDelay = dVal; val.carti_roiDelay = dVal;
[val.cartr_roiDelay.x,val.carti_roiDelay.x] = pol2cart(wrapToPi(valPhase-repmat(phaseShift_roi,size(valPhase,1),1)),valMag./repmat(magScale,size(valMag,1),1));

info = fields(tr)';
for infoInd = 1:length(info)
    dTr2{infoInd} = tr.(info{infoInd});
    dTr2{infoInd}.info = info{infoInd};
    dTr2{infoInd}.infoDim = 'vox';
    
    dTe2{infoInd} = te.(info{infoInd});
    dTe2{infoInd}.info = info{infoInd};
    dTe2{infoInd}.infoDim = 'vox';
    
    dVal2{infoInd} = val.(info{infoInd});
    dVal2{infoInd}.info = info{infoInd};
    dVal2{infoInd}.infoDim = 'vox';
end
clear tr te val


function [dTr,dVal,dTe] = trainNtest(dTr,dVal,dTe,p,ALflag)
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
    
    xVal = dVal{i}.x;
%     xVal = xVal-repmat(xShift,[size(xVal,1) 1]);
%     xVal = xVal./repmat(xScale,[size(xVal,1) 1]);
    
    xTe = dTe{i}.x;
%     xTe = xTe-repmat(xShift,[size(xTe,1) 1]);
%     xTe = xTe./repmat(xScale,[size(xTe,1) 1]);
    
    %Train
    dTr{i}.svmStruct = svmtrain(dTr{i}.y,xTr,['-s 0 -t 0 -c ' num2str(p.C) ' -q']);
    [dTr{i}.pred,acc,dTr{i}.yhat] = svmpredict(dTr{i}.y,xTr,dTr{i}.svmStruct,'-q');
    dTr{i}.acc = acc(1);
    
%     %
%     [dTr{i}.w,dTr{i}.a] = getAandW(dTr{i}.svmStruct,dTr{i}.x);
%     tmpxTr = mean(dTr{i}.x(dTr{i}.y==1,:),1)-mean(dTr{i}.x(dTr{i}.y==2,:),1);
%     tmpxVal = mean(dVal{i}.x(dVal{i}.y==1,:),1)-mean(dVal{i}.x(dVal{i}.y==2,:),1);
%     tmpxTe = mean(dTe{i}.x(dTe{i}.y==1,:),1)-mean(dTe{i}.x(dTe{i}.y==2,:),1);
%     [dTr{i}.tmp.RHO_tr2te,dTr{i}.tmp.PVAL_tr2te] = corr(tmpxTr',((tmpxTe+tmpxVal)./2)');
%     [dTr{i}.tmp.RHO_tr2w,dTr{i}.tmp.PVAL_tr2w] = corr(tmpxTr',dTr{i}.w');
%     [dTr{i}.tmp.RHO_tr2a,dTr{i}.tmp.PVAL_tr2a] = corr(tmpxTr',dTr{i}.a');
%     [dTr{i}.tmp.RHO_te2w,dTr{i}.tmp.PVAL_te2w] = corr(((tmpxTe+tmpxVal)./2)',dTr{i}.w');
%     [dTr{i}.tmp.RHO_te2a,dTr{i}.tmp.PVAL_te2a] = corr(((tmpxTe+tmpxVal)./2)',dTr{i}.a');
%     dTr{i}.tmp.N = length(tmpxTr);
    
    %Test validation data (for anti-learning detection)
    [dVal{i}.pred,acc,dVal{i}.yhat] = svmpredict(dVal{i}.y,xVal,dTr{i}.svmStruct,'-q');
    dVal{i}.acc = acc(1);
    if ALflag
        dTr{i}.ALrev = sign(mean(dVal{i}.yhat(dVal{i}.y==1))-mean(dVal{i}.yhat(dVal{i}.y==2)));
%         dTr{i}.ALrev = ((mean(dVal{i}.yhat(dVal{i}.y==1))-mean(dVal{i}.yhat(dVal{i}.y==2))>=0)-0.5)*2;
    else
        dTr{i}.ALrev = 0;
    end
    
    %Test test data (for cross-validation)
    [dTe{i}.pred,acc,dTe{i}.yhat] = svmpredict(dTe{i}.y,xTe,dTr{i}.svmStruct,'-q');
    dTe{i}.acc = acc(1);
    
    %Revert anti-learning to learning
    if ALflag
        dTe{i}.yhat = dTe{i}.yhat*dTr{i}.ALrev;
        dTe{i}.pred = ((dTe{i}.pred-1.5)*dTr{i}.ALrev)+1.5;
        dTe{i}.acc = ((dTe{i}.acc-50)*dTr{i}.ALrev)+50;
        
        dVal{i}.yhat = dVal{i}.yhat*dTr{i}.ALrev;
        dVal{i}.pred = ((dVal{i}.pred-1.5)*dTr{i}.ALrev)+1.5;
        dVal{i}.acc = ((dVal{i}.acc-50)*dTr{i}.ALrev)+50;
    end
end

function [dTr2,dVal2,dTe2] = voxelSplit(d,dTr,dVal,dTe,trInd,valInd,teInd,combInfo,infoStr)
for i  =1:length(dTr{1})
    allInfo{i} = dTr{1}{i}.info;
end
for i  =1:length(dTr)
        dTr2{i} = dTr{i}(ismember(allInfo,combInfo));
        dVal2{i} = dVal{i}(ismember(allInfo,combInfo));
        dTe2{i} = dTe{i}(ismember(allInfo,combInfo));
end
dTr = dTr2; clear dTr2
dVal = dVal2; clear dVal2
dTe = dTe2; clear dTe2

for vox = 1:size(dTr{1}{1}.x,2)
    newTr = nan(size(d.y,1),length(dTr),length(dTr{1}));
    newVal = nan(size(d.y,1),length(dTr),length(dTr{1}));
    newTe = nan(size(d.y,1),length(dTr),length(dTr{1}));
    
    for kSubInd = 1:length(dTr)
        for infoInd = 1:length(dTr{1})
            newTr(trInd{kSubInd},kSubInd,infoInd) = dTr{kSubInd}{infoInd}.x(:,vox);
            newVal(valInd{kSubInd},kSubInd,infoInd) = dVal{kSubInd}{infoInd}.x(:,vox);
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
        
        dVal2{kSubInd}{vox} = dVal{kSubInd}{1};
        dVal2{kSubInd}{vox}.pred = [];
        dVal2{kSubInd}{vox}.yhat = [];
        dVal2{kSubInd}{vox}.acc = [];
        dVal2{kSubInd}{vox}.x = squeeze(newVal(valInd{kSubInd},kSubInd,:));
        dVal2{kSubInd}{vox}.info = infoStr;
        dVal2{kSubInd}{vox}.infoDim = combInfo;
        dVal2{kSubInd}{vox}.subLevel = dVal{kSubInd};
        
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


function [dTr,dVal,dTe,trInd,valInd,teInd] = voxelSplitNcomb(d,p,crossVal,subCrossVal,info,infoComb)
for kSubInd = 1:length(subCrossVal)
    %% Split train and test trials
    [dTr{1,kSubInd},dVal{1,kSubInd},dTe{1,kSubInd},trInd{1,kSubInd},valInd{1,kSubInd},teInd{1,kSubInd},~] = splitTrainValTest(d,crossVal,subCrossVal(kSubInd));
    
    %% SVM on low-level
    % At each voxel, combine info and train on info
    
    for combInd = 1:length(infoComb)
        [dTr_tmp(1,1),dVal_tmp(1,1),dTe_tmp(1,1)] = voxelSplit(d,dTr(1,kSubInd),dVal(1,kSubInd),dTe(1,kSubInd),trInd(1,kSubInd),valInd(1,kSubInd),teInd(1,kSubInd),info(combInd,:),infoComb{combInd});
        [dTr_tmp{1,1},dVal_tmp{1,1},dTe_tmp{1,1}] = trainNtest(dTr_tmp{1,1},dVal_tmp{1,1},dTe_tmp{1,1},p,0);
        dTr_tmp2{combInd}(1,kSubInd) = dTr_tmp;
        dVal_tmp2{combInd}(1,kSubInd) = dVal_tmp;
        dTe_tmp2{combInd}(1,kSubInd) = dTe_tmp;
    end
end
clear dTr_tmp dVal_tmp dTe_tmp dTr dVal dTe
% [dTr_tmp,dVal_tmp,dTe_tmp] = levelUp(d,dTr_tmp,dVal_tmp,dTe_tmp,trInd,valInd,teInd,[],infoComb,0);
for combInd = 1:length(infoComb)
    [dTr_tmp,dVal_tmp,dTe_tmp] = levelUp(d,dTr_tmp2{combInd},dVal_tmp2{combInd},dTe_tmp2{combInd},trInd,valInd,teInd,[],infoComb{combInd},0);
    for kSubInd = 1:length(subCrossVal)
        dTr{1,kSubInd}(combInd) = dTr_tmp{1,kSubInd};
        dVal{1,kSubInd}(combInd) = dVal_tmp{1,kSubInd};
        dTe{1,kSubInd}(combInd) = dTe_tmp{1,kSubInd};
    end
end

function [dTr2,dVal2,dTe2] = catInfo(d,dTr,dVal,dTe,trInd,valInd,teInd,combInfo,infoStr)
if ~isempty(combInfo)
    for i  =1:length(dTr{1})
        allInfo{i} = dTr{1}{i}.info;
    end
    for i  =1:length(dTr)
        dTr2{i} = dTr{i}(ismember(allInfo,combInfo));
        dVal2{i} = dVal{i}(ismember(allInfo,combInfo));
        dTe2{i} = dTe{i}(ismember(allInfo,combInfo));
    end
    dTr = dTr2; clear dTr2
    dVal = dVal2; clear dVal2
    dTe = dTe2; clear dTe2
end

newTr = nan(size(d.y,1),length(dTr),length(dTr{1})*size(d.x,2));
newVal = nan(size(d.y,1),length(dTr),length(dTr{1})*size(d.x,2));
newTe = nan(size(d.y,1),length(dTr),length(dTr{1})*size(d.x,2));
for kSubInd = 1:length(dTr)
    for infoInd = 1:length(dTr{1})
        curInd = (infoInd-1)*size(d.x,2)+1:(infoInd)*size(d.x,2);
        newTr(trInd{kSubInd},kSubInd,curInd) = permute(dTr{kSubInd}{infoInd}.x,[1 3 2]);
        newVal(valInd{kSubInd},kSubInd,curInd) = permute(dVal{kSubInd}{infoInd}.x,[1 3 2]);
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
    
    dVal2{kSubInd}{1} = dVal{kSubInd}{1};
    dVal2{kSubInd}{1}.pred = [];
    dVal2{kSubInd}{1}.yhat = [];
    dVal2{kSubInd}{1}.acc = [];
    dVal2{kSubInd}{1}.x = squeeze(newVal(valInd{kSubInd},kSubInd,:));
    dVal2{kSubInd}{1}.info = infoStr;
    dVal2{kSubInd}{1}.infoDim = combInfo;
    dVal2{kSubInd}{1}.subLevel = dTe{kSubInd};

    dTe2{kSubInd}{1} = dTe{kSubInd}{1};
    dTe2{kSubInd}{1}.pred = [];
    dTe2{kSubInd}{1}.yhat = [];
    dTe2{kSubInd}{1}.acc = [];
    dTe2{kSubInd}{1}.x = squeeze(newTe(teInd{kSubInd},kSubInd,:));
    dTe2{kSubInd}{1}.info = infoStr;
    dTe2{kSubInd}{1}.infoDim = combInfo;
    dTe2{kSubInd}{1}.subLevel = dTe{kSubInd};
end


function [dTr2,dVal2,dTe2] = levelUp(d,dTr,dVal,dTe,trInd,valInd,teInd,combInfo,infoStr,keepSubLevel)
if ~exist('keepSubLevel','var')
    keepSubLevel = 1;
end
if ~isempty(combInfo)
    for i  =1:length(dTr{1})
        allInfo{i} = dTr{1}{i}.info;
    end
    for i  =1:length(dTr)
        dTr2{i} = dTr{i}(ismember(allInfo,combInfo));
        dVal2{i} = dVal{i}(ismember(allInfo,combInfo));
        dTe2{i} = dTe{i}(ismember(allInfo,combInfo));
    end
    dTr = dTr2; clear dTr2
    dVal = dVal2; clear dVal2
    dTe = dTe2; clear dTe2
end


newVal = nan(size(d.y,1),length(dTr),length(dTr{1}));
newTe = nan(size(d.y,1),length(dTr),length(dTr{1}));

for kSubInd = 1:length(dTr)
    for infoInd = 1:length(dTr{1})
        newVal(valInd{kSubInd},kSubInd,infoInd) = dVal{kSubInd}{infoInd}.yhat;
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
    dTr2{kSubInd}{1}.x = squeeze(nanmean(newVal(trInd{kSubInd},:,:),2)); % Not actually averaging across subFolds, rather compiling out the nans
    dTr2{kSubInd}{1}.info = infoStr;
    dTr2{kSubInd}{1}.infoDim = combInfo;
    if keepSubLevel
        dTr2{kSubInd}{1}.subLevel = dTr{kSubInd};
    elseif isfield(dTr2{kSubInd}{1},'subLevel')
        dTr2{kSubInd}{1} = rmfield(dTr2{kSubInd}{1},'subLevel');
    end
    
    dVal2{kSubInd}{1} = dVal{kSubInd}{1};
    dVal2{kSubInd}{1}.pred = [];
    dVal2{kSubInd}{1}.yhat = [];
    dVal2{kSubInd}{1}.acc = [];
    dVal2{kSubInd}{1}.x = squeeze(newVal(valInd{kSubInd},kSubInd,:));
    dVal2{kSubInd}{1}.info = infoStr;
    dVal2{kSubInd}{1}.infoDim = combInfo;
    if keepSubLevel
        dVal2{kSubInd}{1}.subLevel = dVal{kSubInd};
    elseif isfield(dVal2{kSubInd}{1},'subLevel')
        dVal2{kSubInd}{1} = rmfield(dVal2{kSubInd}{1},'subLevel');
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
                case {'svmStruct','acc','ALrev'}
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

