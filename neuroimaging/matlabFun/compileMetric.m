function mOut = compileMetric(mOut,mIn)

for anaInd = 1:length(mIn)
    allFields = fields(mIn{anaInd});
    for fieldInd = 1:length(allFields)
        switch allFields{fieldInd}
            case 'info'
            case {'acc','auc','accAbs','dist','distAbs','distT','distTabs','n','nCorrect','yhat1','yhat2'}
                    mOut{anaInd}.(allFields{fieldInd}) = cat(1,mOut{anaInd}.(allFields{fieldInd}),mIn{anaInd}.(allFields{fieldInd}));
            case 'subLevel'
                %deal with different number of sublevels in different
                %sessions (only for the pca analysis)
                if length(mOut{anaInd}.(allFields{fieldInd}))>length(mIn{anaInd}.(allFields{fieldInd}))
                    for i = length(mIn{anaInd}.(allFields{fieldInd}))+1:length(mOut{anaInd}.(allFields{fieldInd}))
                        tmp = mOut{anaInd}.(allFields{fieldInd}){i};
                        allFields2 = fields(tmp);
                        for ii = 1:length(allFields2)
                            switch allFields{ii}
                                case {'auc' 'acc' 'accAbs' 'dist' 'distAbs' 'distT' 'distTabs','yhat1','yhat2'}
                                    tmp.(allFields{ii}) = nan(size(tmp.(allFields{ii})));
                                case {'n' 'nCorrect'}
                                    tmp.(allFields{ii}) = zeros(size(tmp.(allFields{ii})));
                                case 'info'
                                otherwise
                                    error('XX')
                            end
                        end
                        mIn{anaInd}.(allFields{fieldInd}){i} = tmp;
                    end
                elseif length(mOut{anaInd}.(allFields{fieldInd}))<length(mIn{anaInd}.(allFields{fieldInd}))
                    for i = length(mOut{anaInd}.(allFields{fieldInd}))+1:length(mIn{anaInd}.(allFields{fieldInd}))
                        tmp = mIn{anaInd}.(allFields{fieldInd}){i};
                        allFields2 = fields(tmp);
                        for ii = 1:length(allFields2)
                            switch allFields{ii}
                                case {'auc' 'acc' 'accAbs' 'dist' 'distAbs' 'distT' 'distTabs','yhat1','yhat2'}
                                    tmp.(allFields{ii}) = nan(size(tmp.(allFields{ii})));
                                case {'n' 'nCorrect'}
                                    tmp.(allFields{ii}) = zeros(size(tmp.(allFields{ii})));
                                case 'info'
                                otherwise
                                    error('XX')
                            end
                        end
                        mOut{anaInd}.(allFields{fieldInd}){i} = tmp;
                    end
                end
                %go in sub-level
                mOut{anaInd}.(allFields{fieldInd}) = compileMetric(mOut{anaInd}.(allFields{fieldInd}),mIn{anaInd}.(allFields{fieldInd}));
            otherwise
                error('XX');
        end
    end
end

