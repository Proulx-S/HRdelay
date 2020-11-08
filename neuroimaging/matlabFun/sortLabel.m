function rOut = sortLabel(r)

for anaInd = 1:length(r)
    for rep = 1:size(r{anaInd}.y,2)
        yhat1 = r{anaInd}.yhat(r{anaInd}.y(:,rep)==1,:,:,rep);
        yhat2 = r{anaInd}.yhat(r{anaInd}.y(:,rep)==2,:,:,rep);
        rOut{anaInd}.yhat(:,:,:,rep) = cat(1,yhat1,yhat2);
        
        y1 = r{anaInd}.y(r{anaInd}.y(:,rep)==1,rep);
        y2 = r{anaInd}.y(r{anaInd}.y(:,rep)==2,rep);
        rOut{anaInd}.y(:,rep) = cat(1,y1,y2);
    end
    if isfield(r{anaInd},'subLevel') && isfield(r{anaInd}.subLevel{1},'acc')
        rOut{anaInd}.subLevel = sortLabel(r{anaInd}.subLevel);
    end
    if isfield(r{anaInd},'info')
        rOut{anaInd}.info = r{anaInd}.info;
    end
end