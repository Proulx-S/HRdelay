function crossVal = getRandCrossValFolds(d,p)
labelVal = unique(d.label);


%% Sort data to the original structure (structure before label permutation, this has no impact when done before permutation)

%Sort labels
d.ind = (1:length(d.label))'; % To sorte everything back after k-folds assignement
[~,initSortInd] = sort(d.label);
allFields = fields(d);
for curFieldsInd = 1:length(allFields)
    curField = allFields{curFieldsInd};
    if size(d.(curField),1)==p.nObs
        d.(curField) = d.(curField)(initSortInd,:);
    end
end

%Within each label, sort sessions
for labelValInd = 1:length(labelVal)
    curLabel = labelVal(labelValInd);
    curLabelInd = find(d.label==curLabel);
    [a,b] = sort(d.sessionLabel(curLabelInd));
    allFields = fields(d);
    for curFieldsInd = 1:length(allFields)
        curField = allFields{curFieldsInd};
        if size(d.(curField),1)==p.nObs
            d.(curField)(curLabelInd,:) = d.(curField)(curLabelInd(b),:);
        end
    end
end

%Arranged labels in an interleaved fashion
ind1 = 1:length(d.label)/2;
ind2 = ind1+length(d.label)/2;
ind(1:2:length(d.label)) = ind1;
ind(2:2:length(d.label)) = ind2;
allFields = fields(d);
for curFieldsInd = 1:length(allFields)
    curField = allFields{curFieldsInd};
    if size(d.(curField),1)==p.nObs
        d.(curField) = d.(curField)(ind,:);
    end
end
% allStuff2 = [d.label d.sessionLabel d.ind];

%% Assign k-folds (spread session across cross-validation folds, instead of the reverse)
crossVal = zeros(size(d.label));
sessionVal = unique(d.sessionLabel);

sessionInFold = zeros(length(crossVal)/p.k,p.k);
for i = 1:length(sessionVal)
    nInSession(i) = length(find(d.sessionLabel==sessionVal(i)));
end
k = 1;
n = 1;
while any(any(sessionInFold==0))
    %Find session with the most non-assigned obs
    [~,b] = max(nInSession);
    b = find(nInSession(b)==nInSession);
    %Break tighs with the session that has the least obs assigned to that
    %fold
    for curSession = 1:length(sessionVal)
        curN(curSession) = length(find(sessionInFold(:,k)==curSession));
    end
    [~,c] = min(curN(b));
    b = b(c);
    %Randomly break remaining tighs
    b = b(randperm(length(b),1));
    
    sessionInFold(n:n+1,k) = b;
    nInSession(b) = nInSession(b)-2;
        
    if k==p.k
        k = 1;
        n = n+2;
    else
        k = k+1;
    end
end
[d.label d.sessionLabel d.ind];


%Rearange crossVal, two at a time to keep label balanced
for i = 1:2:length(d.sessionLabel)
    for k = 1:size(sessionInFold,2)
        tmp = find(sessionInFold(:,k)==d.sessionLabel(i));
        if ~isempty(tmp)
            crossVal(i:i+1) = k;
            sessionInFold(tmp(1:2),k) = 0;
            break
        end
    end
end
d.crossVal = crossVal;
% [d.label d.sessionLabel d.crossVal d.ind];

% for k = 1:p.k
%     d.sessionLabel(d.crossVal==k)
% end





%% Resort to the original
[~,resortInd] = sort(d.ind);
allFields = fields(d);
for curFieldsInd = 1:length(allFields)
    curField = allFields{curFieldsInd};
    if size(d.(curField),1)==p.nObs
        d.(curField) = d.(curField)(resortInd,:);
    end
end
% allStuff3 = [d.label d.sessionLabel d.ind];
% allStuff==allStuff3

% for k = 1:p.k
%     d.sessionLabel(d.crossVal==k)
% end



%% Shuffle within each label*session intersections to maintain label*session k-fold structure
labelVal = unique(d.label);
sessionVal = unique(d.sessionLabel);
for labelIndVal = 1:length(labelVal)
    for sessionIndVal = 1:length(sessionVal)
        curInd = intersect(find(d.label==labelVal(labelIndVal)),find(d.sessionLabel==sessionVal(sessionIndVal)));
        curInd_shuf = curInd(randperm(length(curInd)));
        d.crossVal(curInd) = d.crossVal(curInd_shuf);
    end
end
crossVal = d.crossVal;
% [d.label(curInd) d.sessionLabel(curInd) d.crossVal(curInd) d.ind(curInd)];

