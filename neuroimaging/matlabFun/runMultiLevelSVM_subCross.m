function [out,d] = runMultiLevelSVM_subCross(dOrig,p)
postProcessType = 'onePair'; % 'none' 'onePair' 'multiClass'

p.getPattern = 1;
p.perm = p.doPerm;

d.x = cat(3,real(dOrig.xData),imag(dOrig.xData));
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

%% Permutate if needed
if p.perm
    error('X')
    d.x = d.x(randperm(size(d.x,1)),:,:);
end

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
    level2{1}.a = zeros(size(d.x,3),size(d.labelPairs,1));
    level2{1}.w = zeros(size(d.x,3),size(d.labelPairs,1));
    level2{2}.a = zeros(size(d.x,2),size(d.labelPairs,1));
    level2{2}.w = zeros(size(d.x,2),size(d.labelPairs,1));
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
    [trPhase,trMag] = cart2pol(dTr.x(:,:,1),dTr.x(:,:,2));
    [tePhase,teMag] = cart2pol(dTe.x(:,:,1),dTe.x(:,:,2));
%     phaseShift1 = circ_mean(trPhase(1:end/2,:),[],1);
%     phaseShift2 = circ_mean(trPhase(end/2+1:end,:),[],1);
%     phaseShift = circ_mean([phaseShift1; phaseShift2],[],1); clear phaseShift1 phaseShift2
%     phaseShift = circ_median(phaseShift,2);
%     trPhase = wrapToPi(trPhase - repmat(phaseShift,size(trPhase,1),size(trPhase,2)));
%     tePhase = wrapToPi(tePhase - repmat(phaseShift,size(tePhase,1),size(tePhase,2)));
    phaseShift = circ_mean([circ_mean(trPhase(dTr.y==1,:),[],1);circ_mean(trPhase(dTr.y==2,:),[],1)],[],1);
    trPhase = wrapToPi(trPhase - repmat(phaseShift,size(trPhase,1),1));
    tePhase = wrapToPi(tePhase - repmat(phaseShift,size(tePhase,1),1));
    magScale = mean(trMag,1);
    trMag = trMag ./ repmat(magScale,size(trMag,1),1);
    teMag = teMag ./ repmat(magScale,size(teMag,1),1);
    [X,Y] = pol2cart(trPhase,trMag); clear trPhase trMag
    dTr.x = cat(3,X,Y);
    [X,Y] = pol2cart(tePhase,teMag); clear tePhase teMag
    dTe.x = cat(3,X,Y);
    
    
    %% Loop over label pairs
    %%%%%
    %%%%%
    %% Train
    %%%%%
    %%%%%
%     dDim1Tr_ord1 = cell(size(d.labelPairs,1),sz(end));% dDim1Te = cell(size(d1Tr.x,3),1);
%     dDim2Tr = cell(size(d.labelPairs,1),1);
    curInd1 = [];
    curInd2 = [];
    for pairInd = 1:size(d.labelPairs,1)      
        
        %% Hierarchical
        % Prepare data
        curInd1 = dTr.y == d.labelPairs(pairInd,1); curInd2 = dTr.y == d.labelPairs(pairInd,2);
        curInd = curInd1|curInd2;
        d1Tr = dTr;
        d1Tr.y = cat(1,ones(length(find(curInd1)),1),ones(length(find(curInd2)),1).*2);
        d1Tr.x = dTr.x(curInd,:,:);
        d1Tr.crossVal = dTr.crossVal(curInd);
        d1Tr.labelPair = d.labelPairs(pairInd,:);
        
        % Run SVM on freq then chan
        [d1Dim1Tr_ord1(pairInd,:),d1Dim2Tr_ord1(pairInd)] = train2Levels(d1Tr,p);
        
        % Run SVM on chan then freq
        d1Tr.x = permute(d1Tr.x,[1 3 2]);
%         d1Tr.x = cat(1,d1Tr.x(end/2+1:end,:,:),d1Tr.x(1:end/2,:,:));
        [d1Dim1Tr_ord2(pairInd,:),d1Dim2Tr_ord2(pairInd)] = train2Levels(d1Tr,p);
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
        curInd1{pairInd} = dTe.y == d.labelPairs(pairInd,1); curInd2{pairInd} = dTe.y == d.labelPairs(pairInd,2);
        curInd = curInd1{pairInd}|curInd2{pairInd};
        d1Te = dTe;
        d1Te.y = cat(1,ones(length(find(curInd1{pairInd})),1),ones(length(find(curInd2{pairInd})),1).*2);
        d1Te.x = dTe.x(curInd,:,:);
        d1Te.labelPair = d.labelPairs(pairInd,:);
        
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
    
    %% Compile fold results
    for pairInd = 1:size(d.labelPairs,1)
        %Predictions
        yhat(teInd(curInd1{pairInd}|curInd2{pairInd}),d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = d2Te{pairInd}.yhat;
        pred(teInd(curInd1{pairInd}|curInd2{pairInd}),d.labelPairs(pairInd,1),d.labelPairs(pairInd,2)) = d2Te{pairInd}.pred;
        acc(d.labelPairs(pairInd,1),d.labelPairs(pairInd,2),kInd) = d2Te{pairInd}.acc; % labelPair x testLabel x trainLabel x k
        
        
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
            
            level3{1}.a(:,pairInd) = level3{1}.a(:,pairInd) + d2Tr{pairInd}.a';
            level3{1}.w(:,pairInd) = level3{1}.w(:,pairInd) + d2Tr{pairInd}.w';
            level3{1}.info = '(freqThenChan; chanThenFreq) x pair';
            level2{1}.a(:,pairInd) = level2{1}.a(:,pairInd) + d1Dim2Tr_ord1{pairInd}.a';
            level2{1}.w(:,pairInd) = level2{1}.w(:,pairInd) + d1Dim2Tr_ord1{pairInd}.w';
            level2{1}.info = '(freqThenChan -> chan) x pair';
            level2{2}.a(:,pairInd) = level2{2}.a(:,pairInd) + d1Dim2Tr_ord2{pairInd}.a';
            level2{2}.w(:,pairInd) = level2{2}.w(:,pairInd) + d1Dim2Tr_ord2{pairInd}.w';
            level2{2}.info = '(chanThenFreq -> freq) x pair';
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
end

switch postProcessType
    case 'onePair'
        yhat = yhat(:,d.labelPairs(1,1),d.labelPairs(1,2));
        pred = pred(:,d.labelPairs(1,1),d.labelPairs(1,2));
        acc = squeeze(acc(d.labelPairs(1,1),d.labelPairs(1,2),:));
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
if p.getPattern
    out.aNw.level3 = level3;
    out.aNw.level2 = level2;
end

% sum([yhat(1:end/2)>0 ; yhat(end/2+1:end)<=0])./length(yhat)
% sum(pred==d.y)/length(pred)
% mean(acc)





function [dDim1Tr,dDim2Tr] = train2Levels(d1Tr,p)
sz = size(d1Tr.x);
%% First-level
dDim1Tr = cell(1,sz(end));
dDim2Tr = cell(1,1);
for chanInd = 1:sz(end)
    %Get data
    dDim1Tr{chanInd} = d1Tr;
    dDim1Tr{chanInd}.x = d1Tr.x(:,:,chanInd);
    %Do SVM
    dDim1Tr{chanInd} = trainSVMatLevel(dDim1Tr{chanInd},p);
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

function [dDim1Te,dDim2Te] = test2Levels(d1Te,dDim1Tr,dDim2Tr)
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





function d = trainSVMatLevel(d,p)
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
        dTr.norm.shift = mean(dTr.x,1);
        dTr.x = dTr.x-repmat(dTr.norm.shift,size(dTr.x,1),1);
        dTe.x = dTe.x-repmat(dTr.norm.shift,size(dTe.x,1),1);
        dTr.norm.scale = std(dTr.x,[],1);
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
d.norm.shift = mean(d.x,1);
d.x = d.x-repmat(d.norm.shift,size(d.x,1),1);
d.norm.scale = std(d.x,[],1);
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
%Train model
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
