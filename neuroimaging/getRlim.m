function rLim = getRlim(data)

condList = fields(data);
rLim = [];
for condInd = 1:length(condList)
    for i = 1:numel(data.(condList{condInd}))
        rLim = [rLim max(abs(data.(condList{condInd}){i}))];
    end
end
rLim = max(rLim);
