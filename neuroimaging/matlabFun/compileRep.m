function rOut = compileRep(rOut,rIn)

for i = 1:length(rOut)
    allFields = fields(rOut{i});
    for fieldInd = 1:length(allFields)
        switch allFields{fieldInd}
            case {'x','w','a'}
                rOut{i}.(allFields{fieldInd}) = cat(3,rOut{i}.(allFields{fieldInd}),rIn{i}.(allFields{fieldInd}));
            case {'y','crossVal','pred','yhat'}
                rOut{i}.(allFields{fieldInd}) = cat(2,rOut{i}.(allFields{fieldInd}),rIn{i}.(allFields{fieldInd}));
            case {'labelPairs','getPattern','info','infoDim','structInfo'}
            case 'subLevel'
                rOut{i}.(allFields{fieldInd}) = compileRep(rOut{i}.(allFields{fieldInd}),rIn{i}.(allFields{fieldInd}));
            case 'norm'
                rOut{i}.(allFields{fieldInd}) = cat(1,rOut{i}.(allFields{fieldInd}),rIn{i}.(allFields{fieldInd}));
            case 'svmStruct'
                rOut{i}.(allFields{fieldInd}) = cat(3,rOut{i}.(allFields{fieldInd}),rIn{i}.(allFields{fieldInd}));
            case 'acc'
                rOut{i}.(allFields{fieldInd}) = cat(2,rOut{i}.(allFields{fieldInd}),rIn{i}.(allFields{fieldInd}));
            otherwise
                error('X')
        end
    end
    if ~isfield(rOut{i},'structInfo')
        for fieldInd = 1:length(allFields)
            switch allFields{fieldInd}
                case 'x'
                    rOut{i}.structInfo.(allFields{fieldInd}) = 'trial x feat x rep';
                case {'w','a'}
                    rOut{i}.structInfo.(allFields{fieldInd}) = 'fold x feat x rep';
                case {'y','crossVal','pred','yhat'}
                    rOut{i}.structInfo.(allFields{fieldInd}) = 'trial x rep';
                case 'infoDim'
                    rOut{i}.structInfo.(allFields{fieldInd}) = 'feat x 1';
                case 'subLevel'
                case 'norm'
                    rOut{i}.structInfo.(allFields{fieldInd}) = 'rep x fold';
                case 'svmStruct'
                    rOut{i}.structInfo.(allFields{fieldInd}) = 'fold x subFold x rep';
                case 'acc'
                    rOut{i}.structInfo.(allFields{fieldInd}) = '1 x rep';
                case {'labelPairs','getPattern','info'}
                otherwise
                    error('X')
            end
        end
    end
end
