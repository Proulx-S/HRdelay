function [out,d] = runMultiLevelSVM(dOrig,p)
postProcessType = 'onePair'; % 'none' 'onePair' 'multiClass'
% p.getPattern = 1;
p.perm = p.doPerm;

d.x = dOrig.xData;
d.y = dOrig.label;
d.crossVal = dOrig.crossVal;
d.label = [1 2];
d.labelPairs = nchoosek(d.label,2);
d.getPattern = p.getPattern;
clear dOrig

% % %% Select equal number of trials at each condition (deterministic OR randomly for each repetition)
% % allInd = [];
% % for i = 1:length(dOrig.label)
% %     tmpInd = find(dOrig.y==dOrig.label(i));
% %     %     allInd = [allInd; tmpInd(randperm(length(tmpInd),dOrig.minLabelRep))];
% %     allInd = [allInd; tmpInd(1:dOrig.minLabelRep)];
% % end
% d = rmfield(dOrig,{'minLabelRep','labelRep'});
% d.labelPairs = nchoosek(d.label,2);
% % d.labelPairs = d.labelPairs(1,:);
% p.getPattern = d.getPattern;

% d.x = d.x(allInd,:,:);
% d.y = d.y(allInd,:,:);


% %% Define cross-validation folds
% % d.crossVal = repmat(1:length(d.y)/length(d.label),1,length(d.label))';
% d.crossVal = [];
% crossVal = (1:length(d.y)/length(d.label))';
% for curLabel = 1:length(d.label)
%     d.crossVal = [d.crossVal; crossVal(randperm(length(crossVal)))];
% end

% %% Permutate if needed
% if p.perm
%     error('X')
%     d.x = d.x(randperm(size(d.x,1)),:,:);
% end

%% Prepare for random permutation if needed
% if p.perm
%     permInd = (1:size(d.x,1))';% permInd = repmat(permInd,1,size(d.labelPairs,1));
% %     for pairInd = 1:size(d.labelPairs,1)
%         for kInd = 1:length(unique(d.crossVal))
%             if rand(1,1)>0.5
%                 permInd(d.crossVal==kInd&d.y==d.labelPairs(pairInd,1),pairInd) = find(d.crossVal==kInd&d.y==d.labelPairs(pairInd,2));
%                 permInd(d.crossVal==kInd&d.y==d.labelPairs(pairInd,2),pairInd) = find(d.crossVal==kInd&d.y==d.labelPairs(pairInd,1));
%                 
%                 permInd(d.crossVal==kInd&d.y==d.labelPairs(pairInd,1),pairInd) = find(d.crossVal==kInd&d.y==d.labelPairs(pairInd,2));
%                 permInd(d.crossVal==kInd&d.y==d.labelPairs(pairInd,2),pairInd) = find(d.crossVal==kInd&d.y==d.labelPairs(pairInd,1));
%             end
%         end
% %     end
% end

%% Prepare params and output
yhat = nan(length(d.y),length(d.label),length(d.label)); %trial x testLabel x trainLabel
pred = nan(length(d.y),length(d.label),length(d.label));
acc = nan(length(d.label),length(d.label),length(unique(d.crossVal)));
yhatL2O1 = nan(length(d.y),length(d.label),length(d.label)); %trial x testLabel x trainLabel
predL2O1 = nan(length(d.y),length(d.label),length(d.label));
accL2O1 = nan(length(d.label),length(d.label),length(unique(d.crossVal)));
yhatL2O2 = nan(length(d.y),length(d.label),length(d.label)); %trial x testLabel x trainLabel
predL2O2 = nan(length(d.y),length(d.label),length(d.label));
accL2O2 = nan(length(d.label),length(d.label),length(unique(d.crossVal)));
yhatL1O1 = nan(length(d.y),length(d.label),length(d.label)); %trial x testLabel x trainLabel
predL1O1 = nan(length(d.y),length(d.label),length(d.label));
accL1O1 = nan(length(d.label),length(d.label),length(unique(d.crossVal)));
yhatL1O2 = nan(length(d.y),length(d.label),length(d.label)); %trial x testLabel x trainLabel
predL1O2 = nan(length(d.y),length(d.label),length(d.label));
accL1O2 = nan(length(d.label),length(d.label),length(unique(d.crossVal)));

% acc = nan(size(d.labelPairs,1),length(d.label),length(d.label),length(unique(d.crossVal)));
% a = zeros(length(d.label),length(d.label),sz(end));
if p.getPattern
%     level3{1}.a = zeros(2,length(d.label),length(d.label));
%     level3{1}.w = zeros(2,length(d.label),length(d.label));
%     level2{1}.a = zeros(size(d.x,3),length(d.label),length(d.label));
%     level2{1}.w = zeros(size(d.x,3),length(d.label),length(d.label));
%     level2{2}.a = zeros(size(d.x,2),length(d.label),length(d.label));
%     level2{2}.w = zeros(size(d.x,2),length(d.label),length(d.label));
    
    level3{1}.a = zeros(2,size(d.labelPairs,1));
    level3{1}.w = zeros(2,size(d.labelPairs,1));
    switch p.infoComb
        case {'all1','all1x'}
            level2{1}.a = zeros(size(d.x,2),size(d.labelPairs,1));
            level2{1}.w = zeros(size(d.x,2),size(d.labelPairs,1));
            level2{2}.a = zeros(size(d.x,2),size(d.labelPairs,1));
            level2{2}.w = zeros(size(d.x,2),size(d.labelPairs,1));
        case {'delayCart','delayCart2'}
            level2{1}.a = zeros(2,size(d.labelPairs,1));
            level2{1}.w = zeros(2,size(d.labelPairs,1));
            level2{2}.a = zeros(size(d.x,2),size(d.labelPairs,1));
            level2{2}.w = zeros(size(d.x,2),size(d.labelPairs,1));
        case {'alli','alli_cross'}
            level2{1}.a = zeros(2,size(d.labelPairs,1));
            level2{1}.w = zeros(2,size(d.labelPairs,1));
            level2{2}.a = zeros(size(d.x,2),size(d.labelPairs,1));
            level2{2}.w = zeros(size(d.x,2),size(d.labelPairs,1));
            level1{1}.a = zeros(size(d.x,2),size(d.labelPairs,1));
            level1{1}.w = zeros(size(d.x,2),size(d.labelPairs,1));
            level1{2}.a = zeros(size(d.x,2),size(d.labelPairs,1));
            level1{2}.w = zeros(size(d.x,2),size(d.labelPairs,1));
        otherwise
            error('X')
    end
end

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
    switch p.infoComb
        case {'delayPol_ampPol','delayPol_ampPol_cross'}
            %Normalize in polar coordinates
            trPhase = angle(dTr.x); trMag = abs(dTr.x);
            tePhase = angle(dTe.x); teMag = abs(dTe.x);
            phaseShift = circ_mean([circ_mean(trPhase(dTr.y==1,:),[],1);circ_mean(trPhase(dTr.y==2,:),[],1)],[],1);
            trPhase = wrapToPi(trPhase - repmat(phaseShift,size(trPhase,1),1));
            tePhase = wrapToPi(tePhase - repmat(phaseShift,size(tePhase,1),1));
            magScale = mean(trMag,1);
            trMag = trMag ./ repmat(magScale,size(trMag,1),1);
            teMag = teMag ./ repmat(magScale,size(teMag,1),1);
            %Put back in hierarchical data matrix
            dTr.x = cat(3,trPhase,trMag); clear trPhase trMag
            dTe.x = cat(3,tePhase,teMag); clear tePhase teMag
%         case 'delayCarti'
%             %Normalize in polar coordinates
%             trPhase = angle(dTr.x); trMag = ones(size(dTr.x));
%             tePhase = angle(dTe.x); teMag = ones(size(dTe.x));
%             phaseShift = circ_mean([circ_mean(trPhase(dTr.y==1,:),[],1);circ_mean(trPhase(dTr.y==2,:),[],1)],[],1);
%             trPhase = wrapToPi(trPhase - repmat(phaseShift,size(trPhase,1),1));
%             tePhase = wrapToPi(tePhase - repmat(phaseShift,size(tePhase,1),1));
%             %Put back in hierarchical data matrix
%             [~,dTr.x] = pol2cart(trPhase,trMag); clear trPhase trMag
%             [~,dTe.x] = pol2cart(tePhase,teMag); clear tePhase teMag
        case {'delayCart','delayCart2','delayCart_cross'}
            %Normalize in polar coordinates
            trPhase = angle(dTr.x); trMag = ones(size(dTr.x));
            tePhase = angle(dTe.x); teMag = ones(size(dTe.x));
            phaseShift = circ_mean([circ_mean(trPhase(dTr.y==1,:),[],1);circ_mean(trPhase(dTr.y==2,:),[],1)],[],1);
            trPhase = wrapToPi(trPhase - repmat(phaseShift,size(trPhase,1),1));
            tePhase = wrapToPi(tePhase - repmat(phaseShift,size(tePhase,1),1));
            %Put back in hierarchical data matrix
            [X,Y] = pol2cart(trPhase,trMag); clear trPhase trMag
            dTr.x = cat(3,X,Y);
            [X,Y] = pol2cart(tePhase,teMag); clear tePhase teMag
            dTe.x = cat(3,X,Y);
        case {'Cart','Cart_cross'}
            %Normalize in polar coordinates
            trPhase = angle(dTr.x); trMag = ones(size(dTr.x));
            tePhase = angle(dTe.x); teMag = ones(size(dTe.x));
            phaseShift = circ_mean([circ_mean(trPhase(dTr.y==1,:),[],1);circ_mean(trPhase(dTr.y==2,:),[],1)],[],1);
            trPhase = wrapToPi(trPhase - repmat(phaseShift,size(trPhase,1),1));
            tePhase = wrapToPi(tePhase - repmat(phaseShift,size(tePhase,1),1));
            magScale = mean(trMag,1);
            trMag = trMag ./ repmat(magScale,size(trMag,1),1);
            teMag = teMag ./ repmat(magScale,size(teMag,1),1);
            %Put back in hierarchical data matrix
            [X,Y] = pol2cart(trPhase,trMag); clear trPhase trMag
            dTr.x = cat(3,X,Y); clear X Y
            [X,Y] = pol2cart(tePhase,teMag); clear tePhase teMag
            dTe.x = cat(3,X,Y); clear X Y
        case {'all1','all1x','all2'}
            %Normalize in polar coordinates
            trPhase = angle(dTr.x); trMag = abs(dTr.x);
            tePhase = angle(dTe.x); teMag = abs(dTe.x);
            phaseShift = circ_mean([circ_mean(trPhase(dTr.y==1,:),[],1);circ_mean(trPhase(dTr.y==2,:),[],1)],[],1);
            trPhase = wrapToPi(trPhase - repmat(phaseShift,size(trPhase,1),1));
            tePhase = wrapToPi(tePhase - repmat(phaseShift,size(tePhase,1),1));
            magScale = mean(trMag,1);
            trMag = trMag ./ repmat(magScale,size(trMag,1),1);
            teMag = teMag ./ repmat(magScale,size(teMag,1),1);
            %Put back in hierarchical data matrix
            [X,Y] = pol2cart(trPhase,trMag); clear trPhase trMag
            dTr.x = complex(X,Y); clear X Y
            [X,Y] = pol2cart(tePhase,teMag); clear tePhase teMag
            dTe.x = complex(X,Y); clear X Y
        case {'alli','alli_cross'}
            %Normalize in polar coordinates
            trPhase = angle(dTr.x); trMag = abs(dTr.x);
            tePhase = angle(dTe.x); teMag = abs(dTe.x);
            phaseShift = circ_mean([circ_mean(trPhase(dTr.y==1,:),[],1);circ_mean(trPhase(dTr.y==2,:),[],1)],[],1);
            trPhase = wrapToPi(trPhase - repmat(phaseShift,size(trPhase,1),1));
            tePhase = wrapToPi(tePhase - repmat(phaseShift,size(tePhase,1),1));
            magScale = mean(trMag,1);
            trMag = trMag ./ repmat(magScale,size(trMag,1),1);
            teMag = teMag ./ repmat(magScale,size(teMag,1),1);
            %Put back in hierarchical data matrix
            [X,Y] = pol2cart(trPhase,trMag); clear trPhase trMag
            dTr.x = complex(X,Y); clear X Y
            dTr.x = cat(3,imag(dTr.x),abs(dTr.x));
            [X,Y] = pol2cart(tePhase,teMag); clear tePhase teMag
            dTe.x = complex(X,Y); clear X Y
            dTe.x = cat(3,imag(dTe.x),abs(dTe.x));
        otherwise
            error('X')
    end
    
    %% Run multilevel
    switch p.infoComb
        case {'delayPol_ampPol','delayPol_ampPol_cross','delayCart','delayCarti','delayCart2','delayCart_cross','Cart','Cart_cross','alli','alli_cross'}
            if strcmp(p.infoComb,'delayCart2')
                skipZforDelayCart = 1;
            else
                skipZforDelayCart = 0;
            end
            [d2Tr,d1Dim1Tr_ord1,d1Dim2Tr_ord1,d1Dim1Tr_ord2,d1Dim2Tr_ord2,d2Te,d1Dim1Te_ord1,d1Dim2Te_ord1,d1Dim1Te_ord2,d1Dim2Te_ord2,curInd,curInd1,curInd2] = trainNtest_2levels(d,p,dTr,dTe,skipZforDelayCart);
%             [d2Tr,d1Dim2Tr_ord1,d1Dim2Tr_ord2,d2Te,d1Dim2Te_ord1,d1Dim2Te_ord2,curInd,curInd1,curInd2] = trainNtest_2levels(d,p,dTr,dTe,skipZforDelayCart);
%         case 'alli'
%             [d2Tr,d1Dim2Tr_ord1,d1Dim2Tr_ord2,d2Te,d1Dim2Te_ord1,d1Dim2Te_ord2,curInd,curInd1,curInd2] = trainNtest_2levels(d,p,dTr,dTe,skipZforDelayCart);
        case {'all1','all1x'}
            % delayCart infoThenVox (L1L2)
            dTr_cur = dTr;
            [X,Y] = pol2cart(angle(dTr_cur.x),ones(size(dTr_cur.x)));
            dTr_cur.x = cat(3,X,Y);
            dTr_cur.x = permute(dTr_cur.x,[1 3 2]);
            dTe_cur = dTe;
            [X,Y] = pol2cart(angle(dTe_cur.x),ones(size(dTe_cur.x)));
            dTe_cur.x = cat(3,X,Y);
            dTe_cur.x = permute(dTe_cur.x,[1 3 2]);
            if strcmp(p.infoComb,'all1x'); signOutput = 1; else signOutput = 0; end
            [trL2_dim1,teL2_dim1,curInd,curInd1,curInd2] = trainNtest_L1L2(d,p,dTr_cur,dTe_cur,signOutput);
            
%             w = nan(2,length(trL2_dim1.L1));
%             a = nan(2,length(trL2_dim1.L1));
%             for i = 1:length(trL2_dim1.L1)
%                 w(:,i) = trL2_dim1.L1{i}.w;
%                 a(:,i) = trL2_dim1.L1{i}.a;
%             end
%             figure('WindowStyle','docked');
%             scatter(w(1,:),a(1,:)); hold on
%             scatter(w(2,:),a(2,:)); legend({'real' 'imag'})
%             xlabel('w');ylabel('a')
%             figure('WindowStyle','docked');
%             hist(w',50); legend({'real' 'imag'}); xlabel('w')
%             yLim = get(gca,'ylim'); xlim([-2 2]);
%             figure('WindowStyle','docked');
%             hist(a',50); legend({'real' 'imag'}); xlabel('a')
%             ylim(yLim); xlim([-2 2]);
            
            
            % ampPol (L2)
            trL2_dim2.L2 = cell(size(d.labelPairs,1),1);
            trL2_dim2 = repmat(trL2_dim2,size(d.labelPairs,1),1);
            teL2_dim2 = trL2_dim2;
            for pairInd = 1:size(teL2_dim1,1)
                % Train
                [dTr_cur,~,~,~] = getLabelPairData(dTr,d.labelPairs(pairInd,:));
                dTr_cur.x = abs(dTr_cur.x);
                dTr_cur = trainSVMatLevel(dTr_cur,p);
                % Test
                [dTe_cur,~,~,~] = getLabelPairData(dTe,d.labelPairs(pairInd,:));
                dTe_cur.x = abs(dTe_cur.x);
                dTe_cur = testSVMatLevel(dTe_cur,dTr_cur);
                % Compile
                trL2_dim2(pairInd).L2 = {dTr_cur};
                teL2_dim2(pairInd).L2 = {dTe_cur};
            end
            
            % L3
            trL3.L1 = cell(size(d.labelPairs,1),1);
            trL3 = repmat(trL3,size(d.labelPairs,1),1);
            teL3 = trL3;
            for pairInd = 1:size(trL3,1)
                % Train
                [dTr_cur,~,~,~] = getLabelPairData(dTr,d.labelPairs(pairInd,:));
                dTr_cur.x = cat(2,trL2_dim1(pairInd).L2{1}.yhat,trL2_dim2(pairInd).L2{1}.yhat);
                dTr_cur = trainSVMatLevel(dTr_cur,p);
                % Test
                [dTe_cur,~,~,~] = getLabelPairData(dTe,d.labelPairs(pairInd,:));
                dTe_cur.x = cat(2,teL2_dim1(pairInd).L2{1}.yhat,teL2_dim2(pairInd).L2{1}.yhat);
                dTe_cur = testSVMatLevel(dTe_cur,dTr_cur);
                % Compile
                trL3(pairInd).L1 = {dTr_cur};
                teL3(pairInd).L1 = {dTe_cur};
            end
            
            %% For backward compatibility
            for pairInd = 1:size(trL3,1)
                d2Te{pairInd}.yhat = teL3(pairInd).L1{1}.yhat;
                d2Te{pairInd}.pred = teL3(pairInd).L1{1}.pred;
                d2Te{pairInd}.acc = teL3(pairInd).L1{1}.acc;
                
                d1Dim2Te_ord1{pairInd}.yhat = teL2_dim1(pairInd).L2{1}.yhat;
                d1Dim2Te_ord1{pairInd}.pred = teL2_dim1(pairInd).L2{1}.pred;
                d1Dim2Te_ord1{pairInd}.acc = teL2_dim1(pairInd).L2{1}.acc;
                d1Dim2Te_ord2{pairInd}.yhat = teL2_dim2(pairInd).L2{1}.yhat;
                d1Dim2Te_ord2{pairInd}.pred = teL2_dim2(pairInd).L2{1}.pred;
                d1Dim2Te_ord2{pairInd}.acc = teL2_dim2(pairInd).L2{1}.acc;
            end
        case 'all2'
            % delayCart infoThenVox (L1L2 (but use only L1))
            dTr_cur = dTr;
            [X,Y] = pol2cart(angle(dTr_cur.x),ones(size(dTr_cur.x)));
            dTr_cur.x = cat(3,X,Y);
            dTr_cur.x = permute(dTr_cur.x,[1 3 2]);
            dTe_cur = dTe;
            [X,Y] = pol2cart(angle(dTe_cur.x),ones(size(dTe_cur.x)));
            dTe_cur.x = cat(3,X,Y);
            dTe_cur.x = permute(dTe_cur.x,[1 3 2]);
            [trL2_dim1,teL2_dim1,curInd,curInd1,curInd2] = trainNtest_L1L2(d,p,dTr_cur,dTe_cur);
            
            % delayCart + ampPol (L2)
            %ord1=info_then_vox
            %ord2=vox_then_info
            %%%%%%%%%%%%%%%%%%%%%
            %WILL NOT BE FULLY COMPATIBLE FOR MULTICLASS
            %%%%%%%%%%%%%%%%%%%%%
            if length(d.label)>2
                error('not compatible for multiclass')
            end
            dTr_cur = dTr; dTr_cur.x = [];
            dTe_cur = dTe; dTe_cur.x = [];
            for voxInd = 1:size(dTr.x,2)
                dTr_cur.x(:,:,voxInd) = cat(2,trL2_dim1.L1{voxInd}.yhat,abs(dTr.x(:,voxInd)));
                dTe_cur.x(:,:,voxInd) = cat(2,teL2_dim1.L1{voxInd}.yhat,abs(dTe.x(:,voxInd)));
            end
            [d2Tr,~,d1Dim2Tr_ord1,~,d1Dim2Tr_ord2,d2Te,~,d1Dim2Te_ord1,~,d1Dim2Te_ord2,~,~,~] = trainNtest_2levels(d,p,dTr_cur,dTe_cur);
        otherwise
            error('X')
    end

    %% Compile fold results
    for pairInd = 1:size(d.labelPairs,1)
        %Predictions at level3
        yhat(teInd(curInd1{pairInd}|curInd2{pairInd}),d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = d2Te{pairInd}.yhat;
        pred(teInd(curInd1{pairInd}|curInd2{pairInd}),d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = d2Te{pairInd}.pred;
        acc(d.labelPairs(pairInd,1),d.labelPairs(pairInd,2),kInd) = d2Te{pairInd}.acc; % labelPair x testLabel x trainLabel x k
        %Predictions at level2
        yhatL2O1(teInd(curInd1{pairInd}|curInd2{pairInd}),d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = d1Dim2Te_ord1{pairInd}.yhat;
        predL2O1(teInd(curInd1{pairInd}|curInd2{pairInd}),d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = d1Dim2Te_ord1{pairInd}.pred;
        accL2O1(d.labelPairs(pairInd,1),d.labelPairs(pairInd,2),kInd) = d1Dim2Te_ord1{pairInd}.acc; % labelPair x testLabel x trainLabel x k
        yhatL2O2(teInd(curInd1{pairInd}|curInd2{pairInd}),d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = d1Dim2Te_ord2{pairInd}.yhat;
        predL2O2(teInd(curInd1{pairInd}|curInd2{pairInd}),d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = d1Dim2Te_ord2{pairInd}.pred;
        accL2O2(d.labelPairs(pairInd,1),d.labelPairs(pairInd,2),kInd) = d1Dim2Te_ord2{pairInd}.acc; % labelPair x testLabel x trainLabel x k
        
        if strcmp(p.infoComb,'alli') || strcmp(p.infoComb,'alli_cross')
            %Predictions at level1
            yhatL1O1(teInd(curInd1{pairInd}|curInd2{pairInd}),d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = d1Dim1Te_ord1{pairInd,1}.yhat;
            predL1O1(teInd(curInd1{pairInd}|curInd2{pairInd}),d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = d1Dim1Te_ord1{pairInd,1}.pred;
            accL1O1(d.labelPairs(pairInd,1),d.labelPairs(pairInd,2),kInd) = d1Dim1Te_ord1{pairInd,1}.acc; % labelPair x testLabel x trainLabel x k
            yhatL1O2(teInd(curInd1{pairInd}|curInd2{pairInd}),d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = d1Dim1Te_ord1{pairInd,2}.yhat;
            predL1O2(teInd(curInd1{pairInd}|curInd2{pairInd}),d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = d1Dim1Te_ord1{pairInd,2}.pred;
            accL1O2(d.labelPairs(pairInd,1),d.labelPairs(pairInd,2),kInd) = d1Dim1Te_ord1{pairInd,2}.acc; % labelPair x testLabel x trainLabel x k
        end
        
        %a and w
        if p.getPattern
%             level3{1}.a(:,d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = level3{1}.a(:,d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) + d2Tr{pairInd}.a';
%             level3{1}.w(:,d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = level3{1}.w(:,d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) + d2Tr{pairInd}.w';
%             level3{1}.info = 'freqThenChan; chanThenFreq';
%             level2{1}.a(:,d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = level2{1}.a(:,d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) + d1Dim2Tr_ord1{pairInd}.a';
%             level2{1}.w(:,d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = level2{1}.w(:,d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) + d1Dim2Tr_ord1{pairInd}.w';
%             level2{1}.info(pairInd,:,kInd) = 'freqThenChan -> chan';
%             level2{2}.a(:,d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = level2{2}.a(:,d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) + d1Dim2Tr_ord2{pairInd}.a';
%             level2{2}.w(:,d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = level2{2}.w(:,d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) + d1Dim2Tr_ord2{pairInd}.w';
%             level2{2}.info(pairInd,:,kInd) = 'chanThenFreq -> freq';
            

            switch p.infoComb
                case {'all1','all1x'}
                    %for backward compatibility
                    level3{1}.a(:,pairInd) = level3{1}.a(:,pairInd) + trL3.L1{1}.a';
                    level3{1}.w(:,pairInd) = level3{1}.w(:,pairInd) + trL3.L1{1}.w';
                    level3{1}.info = 'delayCart x ampPol';
                    level2{1}.a(:,pairInd) = level2{1}.a(:,pairInd) + trL2_dim1.L2{1}.a';
                    level2{1}.w(:,pairInd) = level2{1}.w(:,pairInd) + trL2_dim1.L2{1}.w';
                    level2{1}.info = 'delayCart x vox';
                    level2{2}.a(:,pairInd) = level2{2}.a(:,pairInd) + trL2_dim2.L2{1}.a';
                    level2{2}.w(:,pairInd) = level2{2}.w(:,pairInd) + trL2_dim2.L2{1}.w';
                    level2{2}.info = 'ampPol x vox';
                case {'delayCart','delayCart2'}
                    level3{1}.a(:,pairInd) = level3{1}.a(:,pairInd) + d2Tr{pairInd}.a';
                    level3{1}.w(:,pairInd) = level3{1}.w(:,pairInd) + d2Tr{pairInd}.w';
                    level3{1}.info = '(freqThenChan; chanThenFreq) x pair';
                    level2{1}.a(:,pairInd) = level2{1}.a(:,pairInd) + d1Dim2Tr_ord1{pairInd}.a';
                    level2{1}.w(:,pairInd) = level2{1}.w(:,pairInd) + d1Dim2Tr_ord1{pairInd}.w';
                    level2{1}.info = '(freqThenChan -> chan) x pair';
                    level2{2}.a(:,pairInd) = level2{2}.a(:,pairInd) + d1Dim2Tr_ord2{pairInd}.a';
                    level2{2}.w(:,pairInd) = level2{2}.w(:,pairInd) + d1Dim2Tr_ord2{pairInd}.w';
                    level2{2}.info = '(chanThenFreq -> freq) x pair';
                case {'alli','alli_cross'}
                    level3{1}.a(:,pairInd) = level3{1}.a(:,pairInd) + d2Tr{pairInd}.a';
                    level3{1}.w(:,pairInd) = level3{1}.w(:,pairInd) + d2Tr{pairInd}.w';
                    level3{1}.info = '(freqThenChan; chanThenFreq) x pair';
                    level2{1}.a(:,pairInd) = level2{1}.a(:,pairInd) + d1Dim2Tr_ord1{pairInd}.a';
                    level2{1}.w(:,pairInd) = level2{1}.w(:,pairInd) + d1Dim2Tr_ord1{pairInd}.w';
                    level2{1}.info = '(freqThenChan -> chan) x pair';
                    level2{2}.a(:,pairInd) = level2{2}.a(:,pairInd) + d1Dim2Tr_ord2{pairInd}.a';
                    level2{2}.w(:,pairInd) = level2{2}.w(:,pairInd) + d1Dim2Tr_ord2{pairInd}.w';
                    level2{2}.info = '(chanThenFreq -> freq) x pair';
                    level1{1}.a(:,pairInd) = level1{1}.a(:,pairInd) +d1Dim1Tr_ord1{pairInd,1}.a';
                    level1{1}.w(:,pairInd) = level1{1}.w(:,pairInd) +d1Dim1Tr_ord1{pairInd,1}.w';
                    level1{1}.info = '';
                    level1{2}.a(:,pairInd) = level1{2}.a(:,pairInd) +d1Dim1Tr_ord1{pairInd,2}.a';
                    level1{2}.w(:,pairInd) = level1{2}.w(:,pairInd) +d1Dim1Tr_ord1{pairInd,2}.w';
                    level1{2}.info = '';
                otherwise
                    error('double-check that')
            end
        end
        
        
        
%         if d.getPattern
%             %Structure
%             a(d.labelPairs(pairInd,1),d.labelPairs(pairInd,2),:) = a(d.labelPairs(pairInd,1),d.labelPairs(pairInd,2),:) + permute(d2Tr{pairInd}.a,[1 3 2]);
%             w(d.labelPairs(pairInd,1),d.labelPairs(pairInd,2),:) = w(d.labelPairs(pairInd,1),d.labelPairs(pairInd,2),:) + permute(d2Tr{pairInd}.w,[1 3 2]);
%         end
    end
    toc
end
if p.getPattern
    level3{1}.a = level3{1}.a./length(unique(d.crossVal));
    level3{1}.w = level3{1}.w./length(unique(d.crossVal));
    level2{1}.a = level2{1}.a./length(unique(d.crossVal));
    level2{1}.w = level2{1}.w./length(unique(d.crossVal));
    level2{2}.a = level2{2}.a./length(unique(d.crossVal));
    level2{2}.w = level2{2}.w./length(unique(d.crossVal));
    if strcmp(p.infoComb,'alli') || strcmp(p.infoComb,'alli_cross')
        level1{1}.a = level1{1}.a./length(unique(d.crossVal));
        level1{1}.w = level1{1}.w./length(unique(d.crossVal));    
        level1{2}.a = level1{2}.a./length(unique(d.crossVal));
        level1{2}.w = level1{2}.w./length(unique(d.crossVal));    
    end
end

switch postProcessType
    case 'onePair'
        yhat = yhat(:,d.labelPairs(1,1),d.labelPairs(1,2));
        pred = pred(:,d.labelPairs(1,1),d.labelPairs(1,2));
        acc = squeeze(acc(d.labelPairs(1,1),d.labelPairs(1,2),:));
        yhatL2O1 = yhatL2O1(:,d.labelPairs(1,1),d.labelPairs(1,2));
        predL2O1 = predL2O1(:,d.labelPairs(1,1),d.labelPairs(1,2));
        accL2O1 = squeeze(accL2O1(d.labelPairs(1,1),d.labelPairs(1,2),:));
        yhatL2O2 = yhatL2O2(:,d.labelPairs(1,1),d.labelPairs(1,2));
        predL2O2 = predL2O2(:,d.labelPairs(1,1),d.labelPairs(1,2));
        accL2O2 = squeeze(accL2O2(d.labelPairs(1,1),d.labelPairs(1,2),:));
        if strcmp(p.infoComb,'alli') || strcmp(p.infoComb,'alli_cross')
            yhatL1O1 = yhatL1O1(:,d.labelPairs(1,1),d.labelPairs(1,2));
            predL1O1 = predL1O1(:,d.labelPairs(1,1),d.labelPairs(1,2));
            accL1O1 = squeeze(accL1O1(d.labelPairs(1,1),d.labelPairs(1,2),:));
            yhatL1O2 = yhatL1O2(:,d.labelPairs(1,1),d.labelPairs(1,2));
            predL1O2 = predL1O2(:,d.labelPairs(1,1),d.labelPairs(1,2));
            accL1O2 = squeeze(accL1O2(d.labelPairs(1,1),d.labelPairs(1,2),:));
        end
    case 'multiClass'
        %Flip matrix
        for i = 1:size(yhat,1)
            %yhat
            tmp = squeeze(yhat(i,:,:));
            tmp(isnan(tmp)) = 0;
            diagInd = logical(diag(ones(size(tmp,2),1)));
            tmp = tmp+-tmp';
            tmp(diagInd) = nan;
            yhat(i,:,:) = tmp;
            %pred
            tmp = (squeeze(pred(i,:,:))-1.5)*2;
            tmp(isnan(tmp)) = 0;
            tmp = tmp+-tmp';
            tmp(tmp==0) = nan;
            pred(i,:,:) = tmp/2+1.5;
        end
        %acc
        acc = mean(acc,3);
        % tmp = (squeeze(acc(i,:,:))-1.5)*2;% labelPair x testLabel x trainLabel x k
        acc(isnan(acc)) = 0;
        acc = acc+acc';
        acc(diag(true(size(acc,1),1))) = nan;
        
        %Take only lines with data
        yhat2 = nan(length(d.label),length(d.label),size(yhat,1)/length(d.label)); % testLabel x trainLabel x trial
        pred2 = nan(length(d.label),length(d.label),size(yhat,1)/length(d.label)); % testLabel x trainLabel x trial
        for i = 1:length(d.label)
            yhat2(i,:,:) = permute(yhat(d.y==d.label(i),i,:),[3 1 2]);
            pred2(i,:,:) = permute(pred(d.y==d.label(i),i,:),[3 1 2]);
        end
        yhat = yhat2; clear yhat2
        pred = pred2; clear pred2 
    otherwise
        error('X')
end

out.yhat = yhat;
out.pred = pred;
out.acc = acc./100;
out.yhatL2O1 = yhatL2O1;
out.predL2O1 = predL2O1;
out.accL2O1 = accL2O1./100;
out.yhatL2O2 = yhatL2O2;
out.predL2O2 = predL2O2;
out.accL2O2 = accL2O2./100;
if strcmp(p.infoComb,'alli') || strcmp(p.infoComb,'alli_cross')
    out.yhatL1O1 = yhatL1O1;
    out.predL1O1 = predL1O1;
    out.accL1O1 = accL1O1./100;
    out.yhatL1O2 = yhatL1O2;
    out.predL1O2 = predL1O2;
    out.accL1O2 = accL1O2./100;
end
if p.getPattern
    out.aNw.level3 = level3;
    out.aNw.level2 = level2;
    if strcmp(p.infoComb,'alli') || strcmp(p.infoComb,'alli_cross')
        out.aNw.level1 = level1;
    end
end

% sum([yhat(1:end/2)>0 ; yhat(end/2+1:end)<=0])./length(yhat)
% sum(pred==d.y)/length(pred)
% mean(acc)



function [d2Tr,d1Dim1Tr_ord1,d1Dim2Tr_ord1,d1Dim1Tr_ord2,d1Dim2Tr_ord2,d2Te,d1Dim1Te_ord1,d1Dim2Te_ord1,d1Dim1Te_ord2,d1Dim2Te_ord2,curInd,curInd1,curInd2] = trainNtest_2levels(d,p,dTr,dTe,skipZforDelayCart)
if ~exist('skipZforDelayCart','var')
    skipZforDelayCart = 0;
end
%% Loop over label pairs
%%%%%
%%%%%
%% Train
%%%%%
%%%%%
%     dDim1Tr_ord1 = cell(size(d.labelPairs,1),sz(end));% dDim1Te = cell(size(d1Tr.x,3),1);
%     dDim2Tr = cell(size(d.labelPairs,1),1);
for pairInd = 1:size(d.labelPairs,1)
    
    %% Hierarchical
    % Prepare data
    [d1Tr,~,~,~] = getLabelPairData(dTr,d.labelPairs(pairInd,:));
    
    % Run SVM on freq then chan
    [d1Dim1Tr_ord1(pairInd,:),d1Dim2Tr_ord1(pairInd)] = train2Levels(d1Tr,p);
    
    % Run SVM on chan then freq
    d1Tr.x = permute(d1Tr.x,[1 3 2]);
    %         d1Tr.x = cat(1,d1Tr.x(end/2+1:end,:,:),d1Tr.x(1:end/2,:,:));
    [d1Dim1Tr_ord2(pairInd,:),d1Dim2Tr_ord2(pairInd)] = train2Levels(d1Tr,p,[],skipZforDelayCart);
    %         figure('WindowStyle','docked');
    %         hist(d1Dim2Tr_ord2{1}.a); xlabel('a')
    %         figure('WindowStyle','docked');
    %         hist(d1Dim2Tr_ord2{1}.w); xlabel('w')
    
    %% Supra-Hierarchical
    d2Tr(pairInd) = {d1Tr};
    d2Tr{pairInd}.x = cat(2,d1Dim2Tr_ord1{pairInd}.yhat,d1Dim2Tr_ord2{pairInd}.yhat);
    d2Tr{pairInd} = trainSVMatLevel(d2Tr{pairInd},p);
    
    %         %% Compile a and w
    %         if p.getPattern
    %             %From high to low levels
    %             level3{1}.a(pairInd,:,kInd) = d2Tr{pairInd}.a;
    %             level3{1}.w(pairInd,:,kInd) = d2Tr{pairInd}.w;
    %             level3{1}.info = 'freqThenChan; chanThenFreq';
    %             level2{1}.a(pairInd,:,kInd) = d1Dim2Tr_ord1{pairInd}.a;
    %             level2{1}.w(pairInd,:,kInd) = d1Dim2Tr_ord1{pairInd}.w;
    %             level2{1}.info(pairInd,:,kInd) = 'freqThenChan -> chan';
    %             level2{2}.a(pairInd,:,kInd) = d1Dim2Tr_ord2{pairInd}.a;
    %             level2{2}.w(pairInd,:,kInd) = d1Dim2Tr_ord2{pairInd}.w;
    %             level2{2}.info(pairInd,:,kInd) = 'chanThenFreq -> freq';
    %         end
end

%%%%%
%%%%%
%% Test
%%%%%
%%%%%
%     display('Testing all pairs')
%     dDim1Te = cell(size(d.labelPairs,1),sz(end));% dDim1Te = cell(size(d1Tr.x,3),1);
%     dDim2Te = cell(size(d.labelPairs,1),1);
curInd1 = cell(1,size(d.labelPairs,1));
curInd2 = cell(1,size(d.labelPairs,1));
for pairInd = 1:size(d.labelPairs,1)
    % Prepare data
    [d1Te,curInd{pairInd},curInd1{pairInd},curInd2{pairInd}] = getLabelPairData(dTe,d.labelPairs(pairInd,:));
    
    %% Hierarchical
    % SVM on Freq then Chanels
    [d1Dim1Te_ord1(pairInd,:),d1Dim2Te_ord1(pairInd)] = test2Levels(d1Te,d1Dim1Tr_ord1(pairInd,:),d1Dim2Tr_ord1);
    
    % SVM on Channels then Freq
    d1Te.x = permute(d1Te.x,[1 3 2]);
    [d1Dim1Te_ord2(pairInd,:),d1Dim2Te_ord2(pairInd)] = test2Levels(d1Te,d1Dim1Tr_ord2(pairInd,:),d1Dim2Tr_ord2);
    
    
    %% Supra-Hierarchical
    d2Te(pairInd) = {d1Te};
    d2Te{pairInd}.x = cat(2,d1Dim2Te_ord1{pairInd}.yhat,d1Dim2Te_ord2{pairInd}.yhat);
    d2Te{pairInd} = testSVMatLevel(d2Te{pairInd},d2Tr{pairInd});
end

function [tr,te,curInd,curInd1,curInd2] = trainNtest_L1L2(d,p,dTr,dTe,signOutput)
if ~exist('signOutput','var')
    signOutput = 0;
end
%% Loop over label pairs
%%%%%
%% Train
%%%%%
tr.L1 = cell(size(d.labelPairs,1),size(dTr.x,3));
tr.L2 = cell(size(d.labelPairs,1),1);
tr = repmat(tr,size(d.labelPairs,1),1);
for pairInd = 1:size(d.labelPairs,1)
    % Prepare data
    [dTr_cur,~,~,~] = getLabelPairData(dTr,d.labelPairs(pairInd,:));
%     curInd1 = dTr.y == d.labelPairs(pairInd,1); curInd2 = dTr.y == d.labelPairs(pairInd,2);
%     curInd = curInd1|curInd2;
%     dTr_cur = dTr;
%     dTr_cur.y = cat(1,ones(length(find(curInd1)),1),ones(length(find(curInd2)),1).*2);
%     dTr_cur.x = dTr.x(curInd,:,:);
%     dTr_cur.crossVal = dTr.crossVal(curInd);
%     dTr_cur.labelPair = d.labelPairs(pairInd,:);
    % Train two levels
    [trL1,trL2] = train2Levels(dTr_cur,p,signOutput);
    tr(pairInd).L1 = trL1;
    tr(pairInd).L2 = trL2;
    
    %         %% Compile a and w
    %         if p.getPattern
    %             %From high to low levels
    %             level3{1}.a(pairInd,:,kInd) = d2Tr{pairInd}.a;
    %             level3{1}.w(pairInd,:,kInd) = d2Tr{pairInd}.w;
    %             level3{1}.info = 'freqThenChan; chanThenFreq';
    %             level2{1}.a(pairInd,:,kInd) = d1Dim2Tr_ord1{pairInd}.a;
    %             level2{1}.w(pairInd,:,kInd) = d1Dim2Tr_ord1{pairInd}.w;
    %             level2{1}.info(pairInd,:,kInd) = 'freqThenChan -> chan';
    %             level2{2}.a(pairInd,:,kInd) = d1Dim2Tr_ord2{pairInd}.a;
    %             level2{2}.w(pairInd,:,kInd) = d1Dim2Tr_ord2{pairInd}.w;
    %             level2{2}.info(pairInd,:,kInd) = 'chanThenFreq -> freq';
    %         end
end

%%%%%
%% Test
%%%%%
te.L1 = cell(size(d.labelPairs,1),size(dTr.x,3));
te.L2 = cell(size(d.labelPairs,1),1);
te = repmat(te,size(d.labelPairs,1),1);
curInd1 = cell(1,size(d.labelPairs,1));
curInd2 = cell(1,size(d.labelPairs,1));
for pairInd = 1:size(d.labelPairs,1)
    % Prepare data
    [dTe_cur,curInd{pairInd},curInd1{pairInd},curInd2{pairInd}] = getLabelPairData(dTe,d.labelPairs(pairInd,:));
%     curInd1{pairInd} = dTe.y == d.labelPairs(pairInd,1); curInd2{pairInd} = dTe.y == d.labelPairs(pairInd,2);
%     curInd = curInd1{pairInd}|curInd2{pairInd};
%     dTe_cur = dTe;
%     dTe_cur.y = cat(1,ones(length(find(curInd1{pairInd})),1),ones(length(find(curInd2{pairInd})),1).*2);
%     dTe_cur.x = dTe.x(curInd,:,:);
%     dTe_cur.labelPair = d.labelPairs(pairInd,:);
    
    % Test 2 levels
    [teL1,teL2] = test2Levels(dTe_cur,trL1(pairInd,:),trL2(pairInd),signOutput);
    te(pairInd).L1 = teL1;
    te(pairInd).L2 = teL2;
end


function [dDim1Tr,dDim2Tr] = train2Levels(d1Tr,p,signOutput,skipZforDelayCart)
if ~exist('signOutput','var')|isempty(signOutput); signOutput = 0; end
if ~exist('skipZforDelayCart','var'); skipZforDelayCart = 0; end
sz = size(d1Tr.x);
%% First-level
dDim1Tr = cell(1,sz(end));
dDim2Tr = cell(1,1);
for chanInd = 1:sz(end)
    %Get data
    dDim1Tr{chanInd} = d1Tr;
    dDim1Tr{chanInd}.x = d1Tr.x(:,:,chanInd);
    %Do SVM
    if signOutput; %Sign the output acording to the sign of the imaginary value
        pCur = p; pCur.getPattern = 1;
        dDim1Tr{chanInd} = trainSVMatLevel(dDim1Tr{chanInd},pCur,skipZforDelayCart);
        if dDim1Tr{chanInd}.w(2)~=0
            dDim1Tr{chanInd}.yhat = sign(dDim1Tr{chanInd}.w(2)).*dDim1Tr{chanInd}.yhat;
        end
    else
        dDim1Tr{chanInd} = trainSVMatLevel(dDim1Tr{chanInd},p,skipZforDelayCart);
    end
end

%% Second-level
%Get data
if dDim1Tr{1}.getPattern
    dDim2Tr{1} = rmfield(dDim1Tr{1},{'norm' 'svmStruct' 'yhat' 'bias' 'w' 'a'});%cell(size(d.x,3),1); dDim1Te = cell(size(d.x,3),1);
else
    dDim2Tr{1} = rmfield(dDim1Tr{1},{'norm' 'svmStruct' 'yhat' 'bias'});%cell(size(d.x,3),1); dDim1Te = cell(size(d.x,3),1);
end
dDim2Tr{1}.x = nan(size(dDim1Tr{1}.yhat,1),sz(end));
for chanInd = 1:length(dDim1Tr)
    dDim2Tr{1}.x(:,chanInd) = dDim1Tr{chanInd}.yhat;
end
%Do SVM
dDim2Tr{1} = trainSVMatLevel(dDim2Tr{1},p);

function [dDim1Te,dDim2Te] = test2Levels(d1Te,dDim1Tr,dDim2Tr,signOutput)
if ~exist('signOutput','var'); signOutput = 0; end
sz = size(d1Te.x);
dDim1Te = cell(1,sz(end));
dDim2Te = cell(1,1);
%% First-level
for chanInd = 1:sz(end)
    %Get data
    dDim1Te{chanInd} = d1Te;
    dDim1Te{chanInd}.x = d1Te.x(:,:,chanInd);
    %Do SVM
    dDim1Te{chanInd} = testSVMatLevel(dDim1Te{chanInd},dDim1Tr{chanInd});
    if signOutput
        if dDim1Tr{chanInd}.w(2)~=0
            dDim1Te{chanInd}.yhat = dDim1Te{chanInd}.yhat.*sign(dDim1Tr{chanInd}.w(2));
        end
    end
end

%% Second-level: SVM on chanel dimension, using output from first-level
%Get data
dDim2Te{1} = rmfield(dDim1Te{chanInd},{'yhat' 'acc'});
dDim2Te{1}.x = nan(size(dDim1Te{chanInd}.yhat,1),sz(end));
for chanInd = 1:length(dDim1Tr)
    dDim2Te{1}.x(:,chanInd) = dDim1Te{chanInd}.yhat;
end
%Do SVM
dDim2Te{1} = testSVMatLevel(dDim2Te{1},dDim2Tr{1});




% function d = trainSVMatLevel(d,p)
% %Normalize
% d.norm.shift = mean(d.x,1);
% d.x = d.x-repmat(d.norm.shift,size(d.x,1),1);
% d.norm.scale = std(d.x,[],1);
% d.x = d.x./repmat(d.norm.scale,size(d.x,1),1);
% %Train model
% d.svmStruct = svmtrain(d.y,d.x,['-s 0 -t 0 -c 1 -q']);
% %Apply model to each data point
% [d.pred, ~, d.yhat] = svmpredict(d.y,d.x,d.svmStruct,'-q');
% d.bias = d.svmStruct.rho;
% %Extract model weigths
% if p.getPattern
%     [d.w,d.a] = getAandW(d.svmStruct,d.x);
% end





function d = trainSVMatLevel(d,p,skipZforDelayCart)
if ~exist('skipZforDelayCart','var'); skipZforDelayCart = 0; end
if p.subCrossVal
    %% Redefine crossVal
    crossVal = unique(d.crossVal);
    for kInd = 1:length(crossVal)
        d.crossVal(d.crossVal==crossVal(kInd)) = kInd;
    end
    clear crossVal
    
    %% Loop over cross-validation folds
    crossVal = unique(d.crossVal);
    d.yhat = nan(size(d.y));
    d.pred = nan(size(d.y));
    d.a = nan(length(crossVal),size(d.x,2));
    d.w = nan(length(crossVal),size(d.x,2));
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
    end
    d.w = mean(d.w,1);
    d.a = mean(d.a,1);
end
%% Do on all data
%Normalize
if ~skipZforDelayCart
    d.norm.shift = mean(d.x,1);
else
    d.norm.shift = zeros(1,size(d.x,2));
end
d.x = d.x-repmat(d.norm.shift,size(d.x,1),1);
if ~skipZforDelayCart
    d.norm.scale = std(d.x,[],1);
else
    d.norm.scale = ones(1,size(d.x,2));
end
d.x = d.x./repmat(d.norm.scale,size(d.x,1),1);
%Train model
d.svmStruct = svmtrain(d.y,d.x,['-s 0 -t 0 -c 1 -q']);
d.bias = d.svmStruct.rho;
if p.subCrossVal
    %Extract model weigths
    if p.getPattern
        [d.wAlt,d.aAlt] = getAandW(d.svmStruct,d.x);
    end
else
    %Apply model to each data point
    [d.pred, ~, d.yhat] = svmpredict(d.y,d.x,d.svmStruct,'-q');
    %Extract model weigths
    if p.getPattern
        [d.w,d.a] = getAandW(d.svmStruct,d.x);
    end
end



function dTe = testSVMatLevel(dTe,dTr)
%Normalize
dTe.x = dTe.x-repmat(dTr.norm.shift,size(dTe.x,1),1);
dTe.x = dTe.x./repmat(dTr.norm.scale,size(dTe.x,1),1);
%Test model
[dTe.pred,acc,dTe.yhat] = svmpredict(dTe.y,dTe.x,dTr.svmStruct,'-q');
dTe.acc = acc(1);


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

function [d,curInd,curInd1,curInd2] = getLabelPairData(d,labelPairs)
curInd1 = d.y == labelPairs(1); curInd2 = d.y == labelPairs(2);
curInd = curInd1|curInd2;
d.y = cat(1,ones(length(find(curInd1)),1),ones(length(find(curInd2)),1).*2);
d.x = d.x(curInd,:,:);
d.crossVal = d.crossVal(curInd);
d.labelPair = labelPairs;
