function [tr,te,trLabel,trSessionLabel,teLabel,teSessionLabel,trLabel1,trLabel2] = splitTrainAndTest(d,curFoldLabel,dataTag)

if ~exist('dataTag','var')
    dataTag = 'xData';
end

trainInd = true(size(d.(dataTag),1),1);
trainInd(d.crossVal==curFoldLabel) = false;
testInd = false(size(d.(dataTag),1),1);
testInd(d.crossVal==curFoldLabel) = true;
tr = d.(dataTag)(trainInd,:);
trLabel = d.label(trainInd);
trSessionLabel = d.sessionLabel(trainInd);
trLabel1 = trLabel==1;
trLabel2 = trLabel==2;
te = d.(dataTag)(testInd,:);
teLabel = d.label(testInd);
teSessionLabel = d.sessionLabel(testInd);
% teLabel1 = teLabel==1;
% teLabel2 = teLabel==2;
