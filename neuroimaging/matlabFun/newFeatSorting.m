function [tr,te,sortInd,mi] = newFeatSorting(tr,te,trLabel1,trLabel2,p)
%Compute MI
mi = computeMI(tr,trLabel1,trLabel2);
[~,sortInd] = sort(mi,'descend');

%Apply sorting
tr = tr(:,sortInd);
te = te(:,sortInd);
