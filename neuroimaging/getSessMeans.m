function sessMeans = getSessMeans(data)

condList = fields(data);
sessMeans = nan([length(condList)+1 size(data.ori1)]);
for sessInd = 1:numel(data.ori1)
    for condInd = 1:length(condList)
        %average runs
        sessMeans(condInd,sessInd) = mean(data.(condList{condInd}){sessInd},1);
    end
    sessMeans(length(condList)+1,sessInd) = mean(sessMeans(ismember(condList,{'ori1' 'ori2'}),sessInd),1);
end
