function [curd,tmpp] = handlePermutation(d,p)

curd = d;
tmpp = p;

%Apply random permutation
for curSession = 1:length(unique(d.sessionLabel))
    shufInd = find(d.sessionLabel==curSession);
    shufInd_shuf = shufInd;
    pairedInd = [randperm(length(shufInd)/2) randperm(length(shufInd)/2)]';
    for paireInd = 1:length(shufInd)/2
        tmpShuf = shufInd(pairedInd==paireInd);
        shufInd_shuf(pairedInd==paireInd,:) = tmpShuf(randperm(length(tmpShuf)));
    end
    %             [d.label(shufInd_shuf) d.label(shufInd)]
    %             length(find(d.label(shufInd_shuf)==1))==length(find(d.label(shufInd_shuf)==2))
    curd.label(shufInd) = d.label(shufInd_shuf);
end


