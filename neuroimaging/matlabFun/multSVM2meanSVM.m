function r = multSVM2meanSVM(r)

for anaInd = 1:length(r)
    if isfield(r{anaInd},'subLevel')
        r{anaInd}.subLevel = multSVM2meanSVM(r{anaInd}.subLevel);
        
        r{anaInd}.info = [r{anaInd}.info '_mean'];
        dims = size(r{anaInd}.yhat);
        for ii = 1:length(r{anaInd}.subLevel)
            if ii==1
                curYhat = r{anaInd}.subLevel{ii}.yhat;
            else
                curYhat = cat(length(dims)+1,curYhat,r{anaInd}.subLevel{ii}.yhat);
            end
        end
        r{anaInd}.yhat = mean(curYhat,length(dims)+1);
    end
    
end
