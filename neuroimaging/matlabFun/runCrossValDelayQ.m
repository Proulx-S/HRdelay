function [delay,crossValDelay,N,edges,binCent] = runCrossValDelayQ(dOrig,p,verbose)

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
binCent2 = binCent;
binCent2(1:p.avBin) = [];
binCent2 = [circ_mean(binCent(1:p.avBin),[],2) binCent2];
binCent2(end-p.avBin+1:end) = [];
binCent2 = [binCent2 circ_mean(binCent(end-p.avBin+1:end),[],2)];

for kInd = 1:length(unique(d.crossVal))
    [dTr{1,1},dTe{1,1},trInd{1,1},teInd{1,1}] = splitTrainTest(d,crossVal(kInd));
    trDelayDiff(kInd,:) = wrapToPi(angle(mean(dTr{1,1}.x(dTr{1,1}.y==2,:),1)) - angle(mean(dTr{1,1}.x(dTr{1,1}.y==1,:),1)));
    teDelayDiff(kInd,:) = wrapToPi(angle(mean(dTe{1,1}.x(dTe{1,1}.y==2,:),1)) - angle(mean(dTe{1,1}.x(dTe{1,1}.y==1,:),1)));
    
    quant = quantile(trDelayDiff(kInd,:),p.nBin-1);
    edges(kInd,:) = [-pi; quant(:); pi];
    [N(kInd,:),~,bin] = histcounts(trDelayDiff(kInd,:),edges(kInd,:));
    binCent(kInd,:) = edges(kInd,1:end-1)+diff(edges(kInd,:))/2;
    
    for i = 1:length(N(kInd,:))
        binInd(i,:) = bin==i;
        crossValDelay(kInd,i) = circ_mean(teDelayDiff(kInd,logical(binInd(i,:))),[],2);
        delay(kInd,i) = circ_mean(trDelayDiff(kInd,logical(binInd(i,:))),[],2);
    end
end

%Average across cross-validation folds
edges = circ_mean(edges,[],1);
% binCent = circ_mean(binCent,[],1);
binWidth = diff(edges);
binCent = edges(1:end-1)+binWidth/2;
    
delay = circ_mean(delay,[],1);
crossValDelay = circ_mean(crossValDelay,[],1);
N = mean(N,1);

%Plot
if verbose
    figure('WindowStyle','docked');
    yyaxis right
    plot(circ_mean(binCent,[],1),circ_mean(delay,[],1),'-o'); hold on
    plot(circ_mean(binCent,[],1),circ_mean(crossValDelay,[],1),'-*');
    yyaxis left
    bar(binCent,N);
end





function [dTr,dTe,trInd,teInd] = splitTrainTest(d,crossVal)

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
