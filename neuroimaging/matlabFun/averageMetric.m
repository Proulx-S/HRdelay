function metric = averageMetric(metric)

for anaInd = 1:length(metric)
    weight = (metric{anaInd}.n./repmat(sum(metric{anaInd}.n,1),size(metric{anaInd}.n,1),1));
    allFields = fields(metric{anaInd});
    for fieldInd = 1:length(allFields)
        switch allFields{fieldInd}
            case {'info','yhat1','yhat2'}
            case {'auc','acc','accAbs','dist','distAbs','distT','distTabs'}
                if size(metric{anaInd}.(allFields{fieldInd}),2)==size(weight,2)
                    metric{anaInd}.(allFields{fieldInd}) = nansum(metric{anaInd}.(allFields{fieldInd}) .* weight,1);
                else
                    metric{anaInd}.(allFields{fieldInd}) = nansum(metric{anaInd}.(allFields{fieldInd}),1);
                end
            case {'n','nCorrect'}
                metric{anaInd}.(allFields{fieldInd}) = sum(metric{anaInd}.(allFields{fieldInd}),1);
            case 'subLevel'
                metric{anaInd}.(allFields{fieldInd}) = averageMetric(metric{anaInd}.(allFields{fieldInd}));
            otherwise
                error('XX');
        end
    end
end

