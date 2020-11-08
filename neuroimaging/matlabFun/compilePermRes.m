function rComp = compilePermRes(rComp,featOld,rNew,featNew)



[~,ind] = ismember(featNew,featOld);

for i = 1:size(rNew,1)
    rCur = nan(1,length(featOld));
    rCur(ind) = rNew(i,:);
    rComp = [rComp; rCur];
end


% Push everything to top of the colums
for i = 1:length(featOld)
    nanInd = isnan(rComp(:,i));
    valueInd = ~nanInd;
    curCol = [rComp(valueInd,i); rComp(nanInd,i)];
    
    rComp(:,i) = curCol;
end

% Shim nan rows
rComp(all(isnan(rComp),2),:) = [];

